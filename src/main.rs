use rust_htslib::bcf::{Reader, Read, record};
use std::env;
use std::collections::HashMap;
    
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
fn make_pos_to_gt_index(vcf_path : &String) -> HashMap<(String, i64), Vec<String>> {
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

// extract the at from the vcf
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

// get the genotype from the VCF
fn get_vcf_gt(record : &record::Record) -> Vec<String> {
    let mut gts = Vec::new();
    let genotypes = record.genotypes().expect("Error reading genotypes");
    // todo: support more than one sample
    let num_samples = record.header().sample_count();
    if num_samples < 1 {
        panic!("no samples!");
    }
    let sample_gt = genotypes.get(0);
    gts.push(sample_gt.to_string());
    gts
}

fn infer_gt(child_ats : &Vec<String>, parent_ats : &Vec<String>) -> Vec<String> {

    // to do: need to hook in actual alleles osomrsoim
    for (i, child_at) in child_ats.iter().enumerate() {
        let mut par_allele : i32 = -1;
        for (j, parent_at) in parent_ats.iter().enumerate() {
            if parent_at.find(child_at).is_some() {
                par_allele = j;
                break;
            }
        }
    }
    vec!["1/1".to_string()]
}

// build up map of vcf id to genotype
// top level: get from the pangenie index (looking up by coordinate because ids don't match)
//            this is based on the assumption that variants are identical between the two vcfs
// other levels: the id is inferred using the AT fields in the variant and its parent
//
fn resolve_genotypes(full_vcf_path : &String,
                     decon_id_to_at : &HashMap<String, Vec<String>>,
                     pg_pos_to_gt : &HashMap<(String, i64), Vec<String>>,
                     id_to_genotype : &mut HashMap<String, Vec<String>>) {

    loop {
        let mut vcf = Reader::from_path(full_vcf_path).expect("Error opening VCF");
        let mut added_count : u64 = 0;
        for (_i, record_result) in vcf.records().enumerate() {
            let record = record_result.expect("Fail to read record");
            let record_id = &record.id();
            let id_string = String::from_utf8_lossy(record_id);
            let id_string_cpy = (*id_string).to_string();
            if id_to_genotype.contains_key(&id_string_cpy) {
                // added in previous iteration
                continue;
            } 
            added_count += 1;
            let rid = record.rid().unwrap();
            let chrom = record.header().rid2name(rid).expect("unable to confird rid to chrom");
            let chrom_string = String::from_utf8_lossy(chrom);
            let chrom_string_cpy = (*chrom_string).to_string();
            let res = pg_pos_to_gt.get(&(chrom_string_cpy, record.pos()));
            if res == None {
                // this site isn't in pangenie, let's see if we can find the parent in
                // the deconstruct vcf
                let mut inferred_gt : Vec<String> = Vec::new();
                match record.info(b"PS").string() {
                    Ok(ps_opt) => {
                        match ps_opt {
                            Some(ps_array) => {
                                let query = id_to_genotype.get(&*String::from_utf8_lossy(ps_array[0]));
                                if query.is_some() {
                                    inferred_gt = infer_gt(&get_vcf_at(&record), &query.unwrap());
                                    println!("scoop {}", String::from_utf8_lossy(ps_array[0]));
                                }
                                
                            }
                            None => ()
                        }
                    }
                    Err(_) => ()                        
                }
                if inferred_gt.len() == 0 {
                    println!("Warning: no parent found for ID {}.  Setting to ./.", id_string_cpy);
                    inferred_gt.push("./.".to_string());
                }
                id_to_genotype.insert(id_string_cpy, inferred_gt);
                println!("Need to look up in decon");
            } else {
                // the site was loaded form pangenie -- just pull the genotype directly
                let gts = &res.unwrap();
                println!("Foind in pangenie: {} => {}", id_string_cpy, gts[0]);
                id_to_genotype.insert(id_string_cpy, gts.to_vec());
            }
        }
        if added_count == 0 {
            break;
        }
    }
}

fn main() -> Result<(), String> {

	 let args: Vec<String> = env::args().collect();
	 let full_vcf_path = &args[1];
    let pg_vcf_path = &args[2];

    // index the full vcf from deconstruct (it must contain all the levels and annotations)
    println!("Indexing deconstruct VCF ATs by ID: {}", full_vcf_path);
    let decon_id_to_at = make_id_to_at_index(full_vcf_path);
    // index the pangenie vcf.  its id's aren't consistent so we use coordinates instead
    println!("Indexing pangenie VCF GTs by position: {}", pg_vcf_path);
    let pg_pos_to_gt = make_pos_to_gt_index(pg_vcf_path);
    
    // this is a map of resolved genotypes
    println!("Resolving genotypes");
    let mut id_to_genotype = HashMap::new();
    resolve_genotypes(full_vcf_path, &decon_id_to_at, &pg_pos_to_gt, &mut id_to_genotype);
    
    Ok(())
}
