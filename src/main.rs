use std::env;
use rust_htslib::bcf::{Reader, Writer, Read, record, Header, Format};
use std::collections::{HashMap, HashSet, VecDeque};
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
}

impl DeconVCFIndex {
    pub fn new() -> DeconVCFIndex {
        DeconVCFIndex {
            id_to_no : HashMap::new(),
            no_to_info : HashMap::new(),
        }
    }

    // add record information to the index
    // called in first pass, so no nesting information is added
    pub fn add(&mut self, id : &str, ats : Vec<String>, alt_allele_names : Vec<String>) -> usize {
        let no : usize = self.id_to_no.len();
        self.id_to_no.insert(id.to_string(), no);
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
        let mut parent_allele_to_child_allele = Vec::new();
        {
            let parent_ats : &Vec<String> = &self.no_to_info.get(parent_no).unwrap().ats;
            let child_ats : &Vec<String> = &self.no_to_info.get(no).unwrap().ats;            
            parent_allele_to_child_allele.resize(parent_ats.len(), -1);
            
            // for every child allele in parallel, compute (child_allele, containing_parent_allele) tuple
            let par_to_chi : Vec<(usize, Option<usize>)> = child_ats
                .par_iter()
                .enumerate()
                .map(|(i, child_at)| {
                    let mut par_idx : Option<usize> = None;
                    // try to find it contained in a parent allele
                    for j in 0..parent_ats.len() {
                        let parent_at : &String = &parent_ats[j];
                        if parent_at.find(child_at).is_some() {
                            // if found, then add the link to the parent's child indexes
                            par_idx = Some(j);
                            break;
                        }
                    }
                    if par_idx.is_none() && false {
                        // if nothing was found, repeat on reversed child traversals
                        let child_rev_at = flip_at(child_at);
                        for j in 0..parent_ats.len() {
                            let parent_at : &String = &parent_ats[j];
                            if parent_at.find(&child_rev_at).is_some() {
                                par_idx = Some(j);
                                break;
                            }
                        }
                    }
                    (i, par_idx)
                })
                .collect();
            // use the tuples to set the parent-to-child index array
            for pc in par_to_chi {
                if pc.1.is_some() {
                    parent_allele_to_child_allele[pc.1.unwrap()] = pc.0 as i16;
                }
            }
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
            for (j, (child_no, child_allele)) in lowest_alleles[i].iter().enumerate() {
                if j > 0 {
                    id_tag.push(':'); // colon-separated list of child alleles for given parent allele
                }
                if child_allele > &0 {
                    id_tag.push_str(&self.no_to_info.get(child_no).unwrap().alt_allele_names[*child_allele as usize - 1]);
                } else {
                    // no child allele that is contained in the parent allele -- so just write it
                    panic!("no child allele found for {}", id);
                }
            }
        }
        id_tag
    }        
}

// Index the sites in the deconstruced VCF
// Done in 2 passes in case any children appear before parents (todo is this even possible?)
fn make_decon_vcf_index(vcf_path : &String) -> DeconVCFIndex {
    let mut index = DeconVCFIndex::new();
    {
        // pass 1: index each site in the deconstructed VCF
        eprintln!("Pass 1 on {}: loading sites into memory", vcf_path);
        let mut vcf = Reader::from_path(vcf_path).expect("Error opening VCF");
        for (_i, record_result) in vcf.records().enumerate() {
            let record = record_result.expect("Fail to read record");
            let record_id = &record.id();
            let id_string = String::from_utf8_lossy(record_id);
            index.add(&id_string, get_vcf_at(&record), get_alt_allele_names(&record));
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
fn make_pos_to_gt_index(vcf_path : &String) -> HashMap<(String, i64), Vec<Vec<i32>>> {
    let mut pos_to_gt = HashMap::new();
    let mut vcf = Reader::from_path(vcf_path).expect("Error opening VCF");    
    for (_i, record_result) in vcf.records().enumerate() {
        let record = record_result.expect("Fail to read record");
        let rid = record.rid().unwrap();
        let chrom = record.header().rid2name(rid).expect("unable to confird rid to chrom");
        let chrom_string = String::from_utf8_lossy(chrom);
        let chrom_string_cpy = (*chrom_string).to_string();
        let gts = get_vcf_gt(&record);
        pos_to_gt.insert((chrom_string_cpy, record.pos()), gts);
    }
    pos_to_gt
}

fn write_resolved_vcf(vcf_path : &String,
                      pg_vcf_path : &String,
                      decon_index : &DeconVCFIndex,
                      gt_index : &HashMap<(String, i64), Vec<Vec<i32>>>) {
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
        out_record.push_info_integer(b"LV", &[get_vcf_level(&in_record)]).expect("Could not set LV");
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
        //Vec<Vec<i32>> genotypes = resolve_genotypes(id_string, decon_index, gt_index);

        let mut gt_vec : Vec<record::GenotypeAllele> = Vec::new();
        for fam in 0..num_samples {
            gt_vec.push(record::GenotypeAllele::Unphased(0));
            gt_vec.push(record::GenotypeAllele::Unphased(0));
        }
        /*
        let found_gt = id_to_genotype.get(&*String::from_utf8_lossy(&out_record.id())).expect("ID not found in map");
        for sample_gt in found_gt {
            for gt_allele in sample_gt {
                if *gt_allele == -1 {
                    gt_vec.push(record::GenotypeAllele::UnphasedMissing);
                } else {
                    gt_vec.push(record::GenotypeAllele::Unphased(*gt_allele));
                }
            }
        }
*/
        out_record.push_genotypes(&gt_vec).expect("Could not set GT");

        out_vcf.write(&out_record).unwrap();
    }        
}

// extract level from the field, returning 0 if not present
fn get_vcf_level(record : &record::Record) -> i32 {
    let lv_info = record.info(b"LV");
    let lv_opt = lv_info.integer().expect("Could not parse LV");
    let lv_int = match lv_opt {
        Some(lv_opt) => {
            let lv_array = *lv_opt;
            lv_array[0]
        },
        None => 0
    };
    lv_int
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
//
// deconstructed VCfs aren't going to fit this mold terribly well (big multiallele sites even in leaves) but we do our best
// it may make more sense just to use snarl ids?
fn get_alt_allele_names(record : &record::Record) -> Vec<String> {
    
    let mut allele_names = Vec::new();
    let alleles = &record.alleles();
    let mut counter = 0 as usize;
    for i in 1..alleles.len() {
        let mut allele_name = String::new();
        let alt_allele = &alleles[i];
        
        let stype;
        if alleles[0].len() == 1 && alt_allele.len() == 1 {
            stype = "SNV".to_string();
        } else if alleles[0].len() == 1 && alt_allele.len() > 1 && alleles[0][0] == alt_allele[0] {
            stype = "INS".to_string();
        } else if alleles[0].len() > 1 && alt_allele.len() == 1 && alleles[0][0] == alt_allele[0] {
            stype = "DEL".to_string();
        } else {
            stype = "COMPLEX".to_string();
        }

        let rid = record.rid().unwrap();
        let chrom = record.header().rid2name(rid).expect("unable to confird rid to chrom");
        let chrom_string = String::from_utf8_lossy(chrom);
        allele_name.push_str(&chrom_string);
        allele_name.push('-');
        // api stores 0-based position which we need to convert to 1-based ourselves
        let pos_1 = record.pos() + 1;
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
            let allele_len = cmp::max(alleles[0].len(), alt_allele.len());
            allele_name.push_str(&allele_len.to_string());
        }
        allele_names.push(allele_name);
    }
    allele_names
}

// parse a string like '0/1' into [0,1]
// dots get turned into -1
fn parse_gt(gt_string : &String) -> Vec<i32> {
    let mut gt_toks : Vec<i32> = Vec::new();
    for gt_tok in gt_string.split(|c| c == '|' || c == '/') {
        if gt_tok == "." {
            gt_toks.push(-1);
        } else {
            gt_toks.push(gt_tok.parse::<i32>().unwrap());
        }        
    }
    gt_toks
}

// get the genotype from the VCF
// each genotype is s vector of ints (where . -> -1)
// and an array is returned with one genotype per sample
fn get_vcf_gt(record : &record::Record) -> Vec<Vec<i32>> {
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

// given a genotype in the parent snarl, infer a genotype in the child snarl
// a child allele gets a parent genotype if its traversal is a substring of
// the parent's traversal.  otherwise it's '.'
fn infer_gt(child_ats : &Vec<String>,
            parent_ats : &Vec<String>,
            parent_gts : &Vec<Vec<i32>>) -> Vec<Vec<i32>> {

    // get every unique parent allele from every sample
    let mut pa_set : HashSet<i32> = HashSet::new();
    for pgt in parent_gts {
        for allele in pgt {
            if *allele != -1 {
                pa_set.insert(*allele);
            }
        }
    }

    // get the reverse complement child traversals
    let mut child_rev_ats : Vec<String> = Vec::new();
    for child_at in child_ats {
        child_rev_ats.push(flip_at(child_at));
    }

    // try to find a child allele for every parent allele, using string comparisons on the traversals
    let mut pa_to_ca : HashMap<i32, i32> = HashMap::new(); // todo: merge with set above
    pa_to_ca.insert(-1, -1);
    for pa in pa_set {
        let p_at : &String = &parent_ats[pa as usize];
        for ca in 0..child_ats.len() {
            if p_at.find(&child_ats[ca]).is_some() ||  p_at.find(&child_rev_ats[ca]).is_some() {
                pa_to_ca.insert(pa, ca as i32);
                break;
            }
        }
        if !pa_to_ca.contains_key(&pa) {
            pa_to_ca.insert(pa, -1);
        }
    }

    // now that we have the allele map, fill in genotypes for the child
    let mut child_gts : Vec<Vec<i32>> = Vec::new();
    for parent_gt in parent_gts {
        let mut child_gt : Vec<i32> = Vec::new();
        for pa in parent_gt {
            child_gt.push(*pa_to_ca.get(pa).unwrap());
        }
        child_gts.push(child_gt);
    }
    child_gts
}

// build up map of vcf id to genotype
// top level: get from the pangenie index (looking up by coordinate because ids don't match)
//            this is based on the assumption that variants are identical between the two vcfs
// other levels: the id is inferred using the AT fields in the variant and its parent
fn resolve_genotypes(full_vcf_path : &String,
                     decon_id_to_at : &HashMap<String, Vec<String>>,
                     pg_pos_to_gt : &HashMap<(String, i64), Vec<Vec<i32>>>,
                     id_to_genotype : &mut HashMap<String, Vec<Vec<i32>>>) {

    // get a null genotype that we'll assign to stuff we don't resolve    
    // todo: clean this
    let mut num_samples : u64 = 1;
    match id_to_genotype.values().next() {
        Some(val) => num_samples = val.len() as u64,
        None => ()
    }
    let mut null_gt : Vec<Vec<i32>> = Vec::new();
    for _i in 0..num_samples {
        null_gt.push(vec![-1,-1]);
    }
    
    let mut added_count_total : u64 = 0;
    let mut loop_count : i32 = 0;
    loop {
        eprintln!("Resolving genotypes [iteration={}]", loop_count);
        let mut vcf = Reader::from_path(full_vcf_path).expect("Error opening VCF");
        let mut added_count : u64 = 0;
        let mut unresolved_ids : Vec<String> = Vec::new();
        let bar = ProgressBar::new(decon_id_to_at.len() as u64);
        let mut log_lines : Vec<String> = Vec::new();
        for (_i, record_result) in vcf.records().enumerate() {
            let record = record_result.expect("Fail to read record");
            let record_id = &record.id();
            let id_string = String::from_utf8_lossy(record_id);
            if id_to_genotype.contains_key(&*id_string) {
                // added in previous iteration
                bar.inc(1);
                continue;
            }
            let id_string_cpy = (*id_string).to_string();
            let rid = record.rid().unwrap();
            let chrom = record.header().rid2name(rid).expect("unable to confird rid to chrom");
            let chrom_string = String::from_utf8_lossy(chrom);
            let chrom_string_cpy = (*chrom_string).to_string();
            let res = pg_pos_to_gt.get(&(chrom_string_cpy, record.pos()));
            let mut gt_found : bool = res != None;
            if gt_found == true {
                // the site was loaded form pangenie -- just pull the genotype directly
                let gts = &res.unwrap();
                id_to_genotype.insert(id_string_cpy, gts.to_vec());
            } else {
                // this site isn't in pangenie, let's see if we can find the parent in
                // the deconstruct vcf
                let mut inferred_gt : Vec<Vec<i32>> = Vec::new();
                let mut parent_found : bool = false;                
                match record.info(b"PS").string() {
                    Ok(ps_opt) => {
                        match ps_opt {
                            Some(ps_array) => {
                                parent_found = true;
                                let parent_gt = id_to_genotype.get(&*String::from_utf8_lossy(ps_array[0]));
                                if parent_gt.is_some() {
                                    let parent_ats = &decon_id_to_at.get(&*String::from_utf8_lossy(ps_array[0])).unwrap();
                                    inferred_gt = infer_gt(&get_vcf_at(&record), parent_ats, &parent_gt.unwrap());
                                    gt_found = true;
                                }                                
                            }
                            None => ()
                        }
                    }
                    Err(_) => ()                        
                }
                if parent_found == false {
                    // no hope of ever genotyping this site if it doesn't have a genotype already and doesn't have a parent record
                    // in the deconstruct vcf.  (this should never happene on the full deconstruct vcf, but could on a subset)
                    log_lines.push(format!("Warning: Top-level (no PS tag) site not present in genotyped vcf (ID {}).  Setting to ./.", id_string_cpy));
                    inferred_gt = null_gt.clone();
                    gt_found = true;
                }
                if gt_found {
                    id_to_genotype.insert(id_string_cpy, inferred_gt);
                }
            }
            if gt_found {
                added_count += 1;
            } else {
                let id_string_cpy = (*id_string).to_string(); // not necessary but rustc can't figure it out
                unresolved_ids.push(id_string_cpy);
            }
            bar.inc(1);
        }
        bar.finish();
        for log_line in log_lines {
            eprintln!("{}", log_line);
        }
        added_count_total += added_count;
        eprintln!("   resolved {} sites", added_count);
        if added_count == 0 || added_count_total == decon_id_to_at.len() as u64 {
            let unresolved_count = unresolved_ids.len();
            for unresolved_id in unresolved_ids {
                eprintln!("Warning: unable to resolve genotype for ID {}. Setting to ./.", unresolved_id);
                id_to_genotype.insert(unresolved_id, null_gt.clone());
            }
            eprintln!("Resolved genotypes for {} / {} sites", decon_id_to_at.len() - unresolved_count, decon_id_to_at.len());
            break;
        }
        loop_count += 1;
    }
}
/*
// scan the deconstructed VCF and print it to stdout, but with the resolved genotypes
fn write_resolved_vcf(full_vcf_path : &String,
                      pg_vcf_path : &String,
                      id_to_genotype : &HashMap<String, Vec<Vec<i32>>>,
                      snarl_forest : &SnarlForest) {

    // we get the samples from the pg vcf
    let pg_vcf = Reader::from_path(pg_vcf_path).expect("Error opening VCF");
    let pg_samples = pg_vcf.header().samples();

    // we write in terms of the original vcf
    // todo: it would be nice to carry over GQ and other fields of interest
    let mut in_vcf = Reader::from_path(full_vcf_path).expect("Error opening VCF");
    let in_header = in_vcf.header();
    let mut out_header = Header::from_template_subset(in_header, &[]).unwrap();
    const REMOVE_TAGS : &'static [&'static str] = &["CONFLICT", "AC", "AF", "NS", "AN", "AT"];
    for info_tag in REMOVE_TAGS {
        out_header.remove_info(info_tag.as_bytes());
    }
    out_header.push_record(b"##INFO=<ID=ID,Number=1,Type=String,Description=\"Colon-separated list of leaf HGSVC-style IDs\"");
    for pg_sample in pg_samples {
        out_header.push_sample(pg_sample);
    }
    let mut out_vcf = Writer::from_stdout(&out_header, true, Format::Vcf).unwrap();
    for (_i, record_result) in in_vcf.records().enumerate() {
        let in_record = record_result.expect("Fail to read record");
        // todo: could be faster to edit in_record in place?  probably doesn't make much difference either way
        let mut out_record = out_vcf.empty_record();
        out_record.set_rid(in_record.rid());
        out_record.set_pos(in_record.pos());
        out_record.set_id(&in_record.id()).expect("Could not set ID");
        out_record.set_alleles(&in_record.alleles()).expect("Could not set alleles");
        out_record.push_info_integer(b"LV", &[get_vcf_level(&in_record)]).expect("Could not set LV");
        match get_vcf_ps(&in_record) {
            Some(ps_string) => out_record.push_info_string(b"PS", &[ps_string.as_bytes()]).expect("Could not set PS"),
            None => (),
        };
        let hprc_id = snarl_forest.get_hprc_id(&String::from_utf8_lossy(&in_record.id()));
        out_record.push_info_string(b"ID", &vec![hprc_id.as_bytes()]).expect("Could not set ID");
        let mut gt_vec : Vec<record::GenotypeAllele> = Vec::new();
        let found_gt = id_to_genotype.get(&*String::from_utf8_lossy(&out_record.id())).expect("ID not found in map");
        for sample_gt in found_gt {
            for gt_allele in sample_gt {
                if *gt_allele == -1 {
                    gt_vec.push(record::GenotypeAllele::UnphasedMissing);
                } else {
                    gt_vec.push(record::GenotypeAllele::Unphased(*gt_allele));
                }
            }
        }
        out_record.push_genotypes(&gt_vec).expect("Could not set GT");
        out_vcf.write(&out_record).unwrap();
    }
}
 */

fn main() -> Result<(), String> {

	 let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: {} <deconstruct vcf> <genotype subset vcf>", &args[0]);
        std::process::exit(1);
    }
	 let full_vcf_path = &args[1];
    let pg_vcf_path = &args[2];
    // todo: cli
    rayon::ThreadPoolBuilder::new().num_threads(8).build_global().unwrap();

    // index the full vcf from deconstruct (it must contain all the levels and annotations)
    let decon_vcf_index = make_decon_vcf_index(full_vcf_path);

    //index the pangenie vcf.  its id's aren't consistent so we use coordinates instead
    eprintln!("Indexing genotyped VCF GTs by position: {}", pg_vcf_path);
    let pg_pos_to_gt = make_pos_to_gt_index(pg_vcf_path);
    
    // this is a map of resolved genotypes
    eprintln!("Writing resolved VCF");
    write_resolved_vcf(full_vcf_path, pg_vcf_path, &decon_vcf_index, &pg_pos_to_gt);
    
    Ok(())
}
