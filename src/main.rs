use clap::{App, Arg};
use rust_htslib::bcf::{Reader, Writer, Read, record, Header, Format};
use std::collections::{HashMap, VecDeque};
use std::cmp;
use indicatif::ProgressBar;
extern crate rayon;
use rayon::prelude::*;

pub struct DeconVCFInfo {
    // link to parent
    parent_no : Option<usize>,
    // key is child node
    // ith value in vector is the index of the allele in the child node that corresponds to the ith allele in this node
    // -1 means it has no allele in the child
    nested_alleles : HashMap<usize, Vec<i16>>,
    // vector of allele traversals
    ats : Vec<String>,
    // verbose hgsvc names, one for each allele
    alt_allele_names : Vec<String>,
}

pub struct DeconVCFIndex {
    // maps ID field (snarl name like >34343>353535) to a number, just to make future references more compact
    id_to_no : HashMap<String, usize>,
    // get the relevant information from a deconstruct VCF record
    no_to_info : HashMap<usize, DeconVCFInfo>,
    // maps position to number (pangenie doesn't seem to keep id for some reason)
    pos_to_no : HashMap<(String, usize), usize>,
    // maps number to pangenie genotype
    no_to_gt : HashMap<usize, Vec<Vec<i16>>>
}

impl DeconVCFIndex {
    pub fn new() -> DeconVCFIndex {
        DeconVCFIndex {
            id_to_no : HashMap::new(),
            no_to_info : HashMap::new(),
            pos_to_no : HashMap::new(),
            no_to_gt : HashMap::new()
        }
    }

    // add record information to the index
    // called in first pass, so no nesting information is added
    pub fn add(&mut self,
               chrom : &str,
               pos : usize,
               id : &str,
               ats : Vec<String>,
               alt_allele_names : Vec<String>) -> usize {
        let no : usize = self.id_to_no.len();
        self.id_to_no.insert(id.to_string(), no);
        self.pos_to_no.insert((chrom.to_string(), pos), no);
        let info = DeconVCFInfo {
            parent_no : None,
            nested_alleles : HashMap::new(),
            ats : ats,
            alt_allele_names : alt_allele_names,
        };
        self.no_to_info.insert(no, info);
        no
    }

    // add nesting information to the index
    // called in second pass, so assume that every site is already in the index
    pub fn add_nesting(&mut self, id : &str, parent_id : &str) -> bool {
        let no = self.id_to_no.get(id).unwrap();
        let parent_no = self.id_to_no.get(parent_id);
        if parent_no.is_none() {
            return false;
        }
        let parent_no = parent_no.unwrap();

        // link from each parent allele to the corresponding contained child allele
        // or -1 if not found
        let parent_allele_to_child_allele : Vec<i16>;
        {
            let parent_ats : &Vec<String> = &self.no_to_info.get(parent_no).unwrap().ats;
            let child_ats : &Vec<String> = &self.no_to_info.get(no).unwrap().ats;
            let child_rev_ats : Vec<String> = child_ats.iter().map(|at| { flip_at(at) }).collect();

            parent_allele_to_child_allele = parent_ats
                .par_iter()
                .map(|parent_at| {
                    for j in 0..child_ats.len() {
                        if parent_at.contains(&child_ats[j]) || parent_at.contains(&child_rev_ats[j]) {
                            return j as i16;
                        }
                    }
                    -1 as i16
                })
                .collect();
        }
            
        //set the parent's allele links to the child alleles
        self.no_to_info
            .get_mut(parent_no)
            .unwrap()
            .nested_alleles
            .insert(*no, parent_allele_to_child_allele);

        //set the child's linke to the parent
        self.no_to_info
            .get_mut(no)
            .unwrap()
            .parent_no = Some(*parent_no);
        true
    }

    // add a genotype
    fn add_genotype(&mut self, chrom : &str, pos : usize, gts : Vec<Vec<i16>>) {
        let no = self.pos_to_no.get(&(chrom.to_string(), pos)).expect(&format!("Could not find ID for site at {} {}", chrom, pos));
        assert_eq!(self.no_to_gt.contains_key(no), false);
        self.no_to_gt.insert(*no, gts);
    }

    // resolve genotypes
    // propagate all genotypes down the allele trees
    // run once after everything is added
    fn resolve_nested_genotypes(&mut self) {
        // do breadth first walk down from all known genotypes
        let mut next_round : Vec<usize> = Vec::new();
        for (no, _) in &self.no_to_gt {
            next_round.push(*no);
        }
        while next_round.len() > 0 {
            let mut next_round_after : Vec<usize> = Vec::new();

            for cur_no in next_round {
                let cur_gts = self.no_to_gt.get(&cur_no).unwrap().clone();
                let nested_alleles = &self.no_to_info.get(&cur_no).unwrap().nested_alleles;

                for (child_no, child_alleles) in nested_alleles {
                    eprintln!("child no {} child alleles {:?}", child_no, child_alleles);
                    let mut child_gts : Vec<Vec<i16>> = Vec::new();
                    for sample_no in 0..cur_gts.len() {
                        let mut child_gt : Vec<i16> = Vec::new();
                        for cur_gt in &cur_gts[sample_no] {
                            let child_al = if *cur_gt >= 0 { child_alleles[*cur_gt as usize] } else { *cur_gt };
                            child_gt.push(child_al);
                        }
                        child_gts.push(child_gt);
                    }
                    self.no_to_gt.insert(*child_no, child_gts);
                    next_round_after.push(*child_no);
                }
            }
            next_round = next_round_after;
        }
        eprintln!("{:?}", self.no_to_gt);
                             
    }

    pub fn get_genotype(&self, id : &str) -> Option<&Vec<Vec<i16>>> {
        let no = self.id_to_no.get(id).unwrap();
        self.no_to_gt.get(no)
    }

    // for a given allele, walk down as far as possible to all child alleles below it
    //
    // todo? goin allele by allele uses way more hash table lookups than we actually need,
    // but doing it that way because it's easiest to think about.  If it's slow in practice
    // either need to adapt datastructure or algorithm...
    fn get_lowest_alleles_for_allele(&self, no : usize, allele : i16) -> Vec<(usize, i16)>{
        // output node/allele pairs
        let mut lowest_alleles : Vec<(usize, i16)> = Vec::new();

        let mut queue : VecDeque<(usize, i16)> = VecDeque::new();
        queue.push_back((no, allele));

        while queue.len() > 0 {
            let (cur_no, cur_allele) = queue.pop_front().unwrap();
            let nested_alleles = &self.no_to_info.get(&cur_no).unwrap().nested_alleles;
            let mut num_pushed = 0;
            if nested_alleles.len() > 0 {
                for (child_no, child_alleles) in nested_alleles {
                    let child_allele = child_alleles[cur_allele as usize];
                    if child_allele >= 0 {
                        queue.push_back((*child_no, child_allele));
                        num_pushed += 1;
                    }
                }
            }
            if num_pushed == 0 {
                // if no alleles found, just push iteself.
                // todo: not sure if this is the desired logic, but the allele is a lowest level allele
                // itself at this point
                lowest_alleles.push((cur_no, cur_allele));
            }
        }
        lowest_alleles        
    }
    
    // return node-allele pair for each allele.  the returned pair corresponds to the lowest node in the
    // tree that has the allele
    fn get_lowest_alleles_for_site(&self, id : &str) -> Vec<Vec<(usize, i16)>> {
        let no = self.id_to_no.get(id).unwrap();
        let info = &self.no_to_info.get(no).unwrap();
        let num_alleles = info.ats.len();
        let mut lowest_alleles : Vec<Vec<(usize, i16)>> = Vec::new();
        for allele in 0..num_alleles {
            lowest_alleles.push(self.get_lowest_alleles_for_allele(*no, allele as i16))
        }
        lowest_alleles
    }

    // make an ID tag
    //
    // quoting Jana:
    //
    // The ID field contains one ID per alternative allele. If an allele contains nested alleles, its ID is composed of a sequence of IDs of the lowest level alleles, separated by “:” (see first ALT allele in example above)
    
    // The IDs itself are unique, for HGSVC we defined them based on the following pattern: 
    // <chrom>-<allele start>-<type>-<REF>-<ALT> (for SNV only) 
    // <chrom>-<allele start>-<type>-<count>-<allele length>  (for other variant types: INS, DEL, COMPLEX)
    pub fn make_id_tag(&self, id : &str) -> String {
        let mut id_tag = String::new();
        let lowest_alleles = self.get_lowest_alleles_for_site(id);
        for i in 1..lowest_alleles.len() { // todo: computed the ref allele children but didn't need to
            if i > 1 {
                id_tag.push(','); // comma-separated list of parent alleles
            }
            eprintln!("lowest alleles {:?}", lowest_alleles[i]);
            let mut ccount : usize = 0;
            for (_j, (child_no, child_allele)) in lowest_alleles[i].iter().enumerate() {
                if child_allele > &0 {
                    if ccount > 0 {
                        id_tag.push(':'); // colon-separated list of child alleles for given parent allele
                    }
                    id_tag.push_str(&self.no_to_info.get(child_no).unwrap().alt_allele_names[*child_allele as usize - 1]);
                    ccount += 1;
                }
            }
        }
        id_tag
    }        
}

// Index the sites in the deconstruced VCF
// Done in 2 passes in case any children appear before parents (todo is this even possible?)
fn make_decon_vcf_index(vcf_path : &str) -> DeconVCFIndex {
    let mut index = DeconVCFIndex::new();
    {
        // pass 1: index each site in the deconstructed VCF
        eprintln!("Pass 1 on {}: loading sites into memory", vcf_path);
        let mut vcf = Reader::from_path(vcf_path).expect("Error opening VCF");
        for (_i, record_result) in vcf.records().enumerate() {
            let record = record_result.expect("Fail to read record");
            let record_id = &record.id();
            let id_string = String::from_utf8_lossy(record_id);
            let rid = record.rid().unwrap();
            let chrom = record.header().rid2name(rid).expect("unable to confird rid to chrom");
            let chrom_string = String::from_utf8_lossy(chrom);
            index.add(&chrom_string, record.pos() as usize, &id_string, get_vcf_at(&record), get_alt_allele_names(&record));
        }        
    }
    {
        // pass 2 add nesting information for every allele
        eprintln!("Pass 2 on {}: loading allele nesting tree", vcf_path);
        let bar = ProgressBar::new(index.id_to_no.len() as u64);        
        let mut vcf = Reader::from_path(vcf_path).expect("Error opening VCF");
        for (_i, record_result) in vcf.records().enumerate() {
            let record = record_result.expect("Fail to read record");
            match get_vcf_ps(&record) {
                Some(ps) => {
                    let record_id = &record.id();
                    let id_string = String::from_utf8_lossy(record_id);
                    index.add_nesting(&id_string, &ps);
                },
                None => ()
            };
            bar.inc(1);
        }
        bar.finish();
    }
    index
}

// Store <CHROM,POS> -> GT for each record in the pangenie VCF
fn make_pos_to_gt_index(vcf_path : &str, index : &mut DeconVCFIndex) {
    let mut vcf = Reader::from_path(vcf_path).expect("Error opening VCF");    
    for (_i, record_result) in vcf.records().enumerate() {
        let record = record_result.expect("Fail to read record");
        let rid = record.rid().unwrap();
        let chrom = record.header().rid2name(rid).expect("unable to confird rid to chrom");
        let chrom_string = String::from_utf8_lossy(chrom);
        let gts = get_vcf_gt(&record);
        index.add_genotype(&chrom_string, record.pos() as usize, gts);
    }
    index.resolve_nested_genotypes();
}

fn write_resolved_vcf(vcf_path : &str,
                      pg_vcf_path : &str,
                      decon_index : &DeconVCFIndex) {
    let mut vcf = Reader::from_path(vcf_path).expect("Error opening VCF");

    // we get the samples from the pg vcf
    let pg_vcf = Reader::from_path(pg_vcf_path).expect("Error opening VCF");
    let pg_samples = pg_vcf.header().samples();

    // we write in terms of the original vcf
    // todo: it would be nice to carry over GQ and other fields of interest
    //       this is done by adding them to the pos_to_gt_index...
    let in_vcf = Reader::from_path(vcf_path).expect("Error opening VCF");
    let in_header = in_vcf.header();
    let mut out_header = Header::from_template_subset(in_header, &[]).unwrap();
    const REMOVE_TAGS : &'static [&'static str] = &["CONFLICT", "AC", "AF", "NS", "AN", "AT"];
    for info_tag in REMOVE_TAGS {
        out_header.remove_info(info_tag.as_bytes());
    }
    out_header.push_record(b"##INFO=<ID=ID,Number=1,Type=String,Description=\"Colon-separated list of leaf HGSVC-style IDs\"");
    let num_samples = pg_samples.len();
    for pg_sample in pg_samples {
        out_header.push_sample(pg_sample);
    }
    let mut null_gt : Vec<record::GenotypeAllele> = Vec::new();
    null_gt.resize(num_samples * 2, record::GenotypeAllele::UnphasedMissing);
    
    let mut out_vcf = Writer::from_stdout(&out_header, true, Format::Vcf).unwrap();

    for (_i, record_result) in vcf.records().enumerate() {
        let in_record = record_result.expect("Fail to read record");
        let record_id = &in_record.id();
        let id_string = String::from_utf8_lossy(record_id);

        // todo: could be faster to edit in_record in place?  probably doesn't make much difference either way
        let mut out_record = out_vcf.empty_record();
        out_record.set_rid(in_record.rid());
        out_record.set_pos(in_record.pos());
        out_record.set_id(&in_record.id()).expect("Could not set ID");
        out_record.set_alleles(&in_record.alleles()).expect("Could not set alleles");
        out_record.push_info_integer(b"LV", &[get_vcf_level(&in_record).into()]).expect("Could not set LV");
        match get_vcf_ps(&in_record) {
            Some(ps_string) => out_record.push_info_string(b"PS", &[ps_string.as_bytes()]).expect("Could not set PS"),
            None => (),
        };

        // get the ID
        // this is a comma-separated list of names for each allele.
        // for non-leaf sites, the name is a colon-separated list of names below
        let id_info_tag = decon_index.make_id_tag(&id_string);
        out_record.push_info_string(b"ID", &vec![id_info_tag.as_bytes()]).expect("Could not set ID");

        // resolve the genotypes
        //Vec<Vec<i16>> genotypes = resolve_genotypes(id_string, decon_index, gt_index);

        let mut gt_vec : Vec<record::GenotypeAllele> = Vec::new();
        
        match decon_index.get_genotype(&String::from_utf8_lossy(&out_record.id())) {
            Some(found_gt) => {
                for sample_gt in found_gt {
                    for gt_allele in sample_gt {
                        if *gt_allele == -1 {
                            gt_vec.push(record::GenotypeAllele::UnphasedMissing);
                        } else {
                            gt_vec.push(record::GenotypeAllele::Unphased(*gt_allele as i32));
                        }
                    }
                }
            },
            None => {
                gt_vec = null_gt.clone();
            }
        };

        out_record.push_genotypes(&gt_vec).expect("Could not set GT");

        out_vcf.write(&out_record).unwrap();
    }        
}

// extract level from the field, returning 0 if not present
fn get_vcf_level(record : &record::Record) -> i16 {
    let lv_info = record.info(b"LV");
    let lv_opt = lv_info.integer().expect("Could not parse LV");
    let lv_int = match lv_opt {
        Some(lv_opt) => {
            let lv_array = *lv_opt;
            lv_array[0]
        },
        None => 0
    };
    lv_int as i16
}

// extract a copy of the ps
fn get_vcf_ps(record : &record::Record) -> Option<String> {
    match record.info(b"PS").string() {
        Ok(ps_opt) => {
            match ps_opt {
                Some(ps_array) => Some(String::from_utf8_lossy(&*ps_array[0]).to_string()),
                None => None,
            }
        },
        Err(_) => None,
    }
}

// extract the AT from the vcf
fn get_vcf_at(record : &record::Record) -> Vec<String> {
    let at_info = record.info(b"AT");
    let mut at_strings = Vec::new();
    let at_res = at_info.string().expect("Could not Parse AT").expect("No AT found");
    let at_array = &*at_res;
    for at_bytes in at_array {
        let s : String = String::from_utf8_lossy(*at_bytes).into_owned();
        at_strings.push(s);
    }
    at_strings
}


// quoting Jana:
//
// The ID field contains one ID per alternative allele. If an allele contains nested alleles, its ID is composed of a sequence of IDs of the lowest level alleles, separated by “:” (see first ALT allele in example above)

// The IDs itself are unique, for HGSVC we defined them based on the following pattern: 
// <chrom>-<allele start>-<type>-<REF>-<ALT> (for SNV only) 
// <chrom>-<allele start>-<type>-<count>-<allele length>  (for other variant types: INS, DEL, COMPLEX)
fn get_alt_allele_names(record : &record::Record) -> Vec<String> {
    
    let mut allele_names = Vec::new();
    let alleles = &record.alleles();
    let mut counter = 0 as usize;
    for i in 1..alleles.len() {
        let mut allele_name = String::new();
        let alt_allele = alleles[i];
        
        let stype;
        let mut slen = 0;
        let mut offset = 0; // if we zip up common prefixes for indel, shift its position to the right a bit
        if alleles[0].len() == 1 && alt_allele.len() == 1 {
            stype = "SNV".to_string();
        } else if alleles[0].len() < alt_allele.len() && &alt_allele[0..alleles[0].len()] == alleles[0] {
            stype = "INS".to_string();
            slen = alt_allele.len() - alleles[0].len();
            offset = alleles[0].len() as i64 - 1;
        } else if alleles[0].len() > alt_allele.len() && &alleles[0][0..alt_allele.len()] == alt_allele {
            stype = "DEL".to_string();
            slen = alleles[0].len() - alt_allele.len();
            offset = alt_allele.len() as i64 - 1;
        } else {
            stype = "COMPLEX".to_string();
            slen = cmp::max(alleles[0].len(), alt_allele.len());
        }

        let rid = record.rid().unwrap();
        let chrom = record.header().rid2name(rid).expect("unable to confird rid to chrom");
        let chrom_string = String::from_utf8_lossy(chrom);
        allele_name.push_str(&chrom_string);
        allele_name.push('-');
        // api stores 0-based position which we need to convert to 1-based ourselves
        let pos_1 = record.pos() + 1 + offset;
        allele_name.push_str(&pos_1.to_string());
        allele_name.push('-');
        allele_name.push_str(&stype);
        if stype == "SNV" {
            allele_name.push('-');
            allele_name.push(alleles[0][0] as char);
            allele_name.push('-');
            allele_name.push(alt_allele[0] as char);
        } else {
            allele_name.push_str(&format!("-{}-", counter));
            counter += 1;
            allele_name.push_str(&slen.to_string());
        }
        allele_names.push(allele_name);
    }
    allele_names
}

// parse a string like '0/1' into [0,1]
// dots get turned into -1
fn parse_gt(gt_string : &String) -> Vec<i16> {
    let mut gt_toks : Vec<i16> = Vec::new();
    for gt_tok in gt_string.split(|c| c == '|' || c == '/') {
        if gt_tok == "." {
            gt_toks.push(-1);
        } else {
            gt_toks.push(gt_tok.parse::<i16>().unwrap());
        }        
    }
    gt_toks
}

// get the genotype from the VCF
// each genotype is s vector of ints (where . -> -1)
// and an array is returned with one genotype per sample
fn get_vcf_gt(record : &record::Record) -> Vec<Vec<i16>> {
    let mut gts = Vec::new();
    let genotypes = record.genotypes().expect("Error reading genotypes");
    // todo: support more than one sample
    let num_samples = record.header().sample_count();
    if num_samples < 1 {
        panic!("no samples!");
    }
    for i in 0..num_samples {
        let sample_gt = genotypes.get(i as usize);
        gts.push(parse_gt(&sample_gt.to_string()));
    }
    gts
}

// turn something like >1>2>3 into <3<2<1
// todo: i think i should be using byte arrays through instead of strings
// which is what the vcf api is doing...
fn flip_at(at : &String) -> String {
    let mut flipped_string : String = String::new();
    flipped_string.reserve_exact(at.len());
    let at_bytes = &at.as_bytes();
    
    let mut i : i64 = at.len() as i64 - 1;
    while i >= 0 {
        // i : last char of step
        // j : first char of step (< or >)
        let mut j : i64 = i - 1;
        while j >= 0 && at_bytes[j as usize] != b'>' && at_bytes[j as usize] != b'<' {
            j -= 1;
        }

        if j >= i || (at_bytes[j as usize] != b'>' && at_bytes[j as usize] != b'<') {
            panic!("Unable to parse AT i={} j={} {}", i, j, at);
        }
        
        if at_bytes[j as usize] == b'>' {
            flipped_string.push('<');
        } else {
            flipped_string.push('>');
        }
        flipped_string.push_str(&at[(j as usize)+1..(i as usize)+1]);
        
        i = j - 1;
    }
    if flipped_string.len() != at.len() {
        panic!("Error flipping {}", at);
    }
    flipped_string
}

fn main() -> Result<(), String> {

    let matches = App::new("resolve-nested-genotypes")
        .version("0.1.0")
        .author("Glenn Hickey <glenn.hickey@gmail.com>")
        .about(concat!(
            "\nPropagte genotypes down to nested alleles using PS and AT tags.\n",
            "The first argument should be a VCF created by vg deconstruct.\n",
            "The second argument should be a genotyped subset of the first. It does not need any tags.\n",
            "Nesting information in LV and PS tags is added back to the output.  An additionl ID tag is\n",
            "added as well, that stores information about the lowest alleles below each allele\n"
        ))
        .arg(
            Arg::with_name("full_vcf_path")
                .value_name("DECONSTRUCT_VCF")
                .help("Full VCF from vg deconstruct")
                .required(true)
                .takes_value(true),
        )
        .arg(
            Arg::with_name("pg_vcf_path")
                .value_name("GENOTYPE_VCF")
                .help("Subset VCF with genotypes")
                .required(true)
                .takes_value(true),
        )
        .arg(
            Arg::with_name("threads")
                .short("t")
                .long("threads")
                .value_name("N")
                .help("Use N threads (default: all available)")
                .takes_value(true),
        )
        .get_matches();
    
	 let full_vcf_path = matches.value_of("full_vcf_path").unwrap();
    let pg_vcf_path = matches.value_of("pg_vcf_path").unwrap();
    match matches.value_of("threads") {
        Some(num_threads) => {
            let num_threads = num_threads.parse::<usize>().expect("Failed to parse threads parameter to int");
            rayon::ThreadPoolBuilder::new().num_threads(num_threads).build_global().unwrap();
        },
        None => ()
    };

    // index the full vcf from deconstruct (it must contain all the levels and annotations)
    let mut decon_vcf_index = make_decon_vcf_index(full_vcf_path);

    // index the pangenie vcf.  its id's aren't consistent so we use coordinates instead
    // also: use the genotypes to resolve every site possible
    eprintln!("Indexing genotyped VCF GTs by position: {}", pg_vcf_path);
    make_pos_to_gt_index(pg_vcf_path, &mut decon_vcf_index);
    
    // this is a map of resolved genotypes
    eprintln!("Writing resolved VCF");
    write_resolved_vcf(full_vcf_path, pg_vcf_path, &decon_vcf_index);
    
    Ok(())
}
