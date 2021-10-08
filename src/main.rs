use std::env;
use rust_htslib::bcf::{Reader, Writer, Read, record, Header, Format};
use std::collections::{HashMap, HashSet};
use indicatif::ProgressBar;
mod snarl_tree;

// Store ID -> AT for each record in the deconstruct VCF
fn make_id_to_at_index(vcf_path : &String) -> HashMap<String, Vec<String>> {
    let mut id_to_at = HashMap::new();
    let mut vcf = Reader::from_path(vcf_path).expect("Error opening VCF");
    for (_i, record_result) in vcf.records().enumerate() {
        let record = record_result.expect("Fail to read record");
        let record_id = &record.id();
        let id_string = String::from_utf8_lossy(record_id);
        let id_string_cpy = (*id_string).to_string();
        let at_strings = get_vcf_at(&record);
        id_to_at.insert(id_string_cpy, at_strings);
    }
    id_to_at
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
                    eprintln!("Warning: Top-level (no PS tag) site not present in genotyped vcf (ID {}).  Setting to ./.", id_string_cpy);
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

// scan the deconstructed VCF and print it to stdout, but with the resolved genotypes
fn write_resolved_vcf(full_vcf_path : &String,
                      pg_vcf_path : &String,
                      id_to_genotype : &HashMap<String, Vec<Vec<i32>>>) {

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

fn main() -> Result<(), String> {

	 let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: {} <deconstruct vcf> <genotype subset vcf>", &args[0]);
        std::process::exit(1);
    }
	 let full_vcf_path = &args[1];
    let pg_vcf_path = &args[2];

    // index the full vcf from deconstruct (it must contain all the levels and annotations)
    eprintln!("Indexing deconstruct VCF ATs by ID: {}", full_vcf_path);
    let decon_id_to_at = make_id_to_at_index(full_vcf_path);

    // index the pangenie vcf.  its id's aren't consistent so we use coordinates instead
    eprintln!("Indexing genotyped VCF GTs by position: {}", pg_vcf_path);
    let pg_pos_to_gt = make_pos_to_gt_index(pg_vcf_path);
    
    // this is a map of resolved genotypes
    let mut id_to_genotype : HashMap<String, Vec<Vec<i32>>> = HashMap::new();
    resolve_genotypes(full_vcf_path, &decon_id_to_at, &pg_pos_to_gt, &mut id_to_genotype);

    // print out the deconstructed VCF but with the genotyped samples
    eprintln!("Writing output");
    write_resolved_vcf(full_vcf_path, pg_vcf_path, &id_to_genotype);
    
    Ok(())
}
