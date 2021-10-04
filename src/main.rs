use vcf::{VCFReader, U8Vec, VCFHeaderFilterAlt, VCFError, VCFRecord};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::BufReader;
use std::env;
use std::collections::HashMap;

// make a ID -> Record in-memory representation of the VCF
fn index_vcf_by_id(vcf_path : &String) -> Result<HashMap<String, VCFRecord>, std::io::Error> {
    let mut mymap = HashMap::new();
    let file = match File::open(vcf_path) {
        Err(why) => panic!("Error opening {}: {}", vcf_path, why),
        Ok(file) => file,
    };
    let decoder = MultiGzDecoder::new(file);
    let buf_reader = BufReader::new(decoder);
    let mut reader = match VCFReader::new(buf_reader) {
        Err(why) => panic!("error parsing vcf {}: {}", vcf_path, why),
        Ok(reader) => reader,
    };

    let mut vcf_record = reader.empty_record();
    loop {
        match reader.next_record(&mut vcf_record) {
            Err(why) => panic!("error iterating vcf {}", why),
            Ok(success) => {
                if success {
                    let id_bytes = vcf_record.id[0].clone();
                    let id_string = String::from_utf8(id_bytes).unwrap();
                    mymap.insert(id_string, vcf_record.clone());                    
                } else {
                    break;
                }
            }
        }
    }
    Ok(mymap)
}

// make a <CHROM,POS> -> Record in-memory representation of the VCF
fn index_vcf_by_pos(vcf_path : &String) -> Result<HashMap<(String, u64), VCFRecord>, std::io::Error> {
    let mut mymap = HashMap::new();
    let file = match File::open(vcf_path) {
        Err(why) => panic!("Error opening {}: {}", vcf_path, why),
        Ok(file) => file,
    };
    let decoder = MultiGzDecoder::new(file);
    let buf_reader = BufReader::new(decoder);
    let mut reader = match VCFReader::new(buf_reader) {
        Err(why) => panic!("error parsing vcf {}: {}", vcf_path, why),
        Ok(reader) => reader,
    };

    let mut vcf_record = reader.empty_record();
    loop {
        match reader.next_record(&mut vcf_record) {
            Err(why) => panic!("error iterating vcf {}", why),
            Ok(success) => {
                if success {
                    let chrom_bytes = vcf_record.chromosome.clone();
                    let chrom_string = String::from_utf8(chrom_bytes).unwrap();
                    mymap.insert((chrom_string, vcf_record.position), vcf_record.clone());                    
                } else {
                    break;
                }
            }
        }
    }
    Ok(mymap)
}

// extract level from the field, returning 0 if not present
fn get_vcf_level(record : &VCFRecord) -> u32 {
    match record.info(b"LV") {
        Some(value) => {
            let lv_bytes = value[0].clone();
            let lv_string = String::from_utf8(lv_bytes).unwrap();
            lv_string.parse().unwrap()
        }
        None => 0,
    }
}

// make a list of record (references) that is sorted decreasing by level
fn sort_by_level<'a>(id_to_record : &'a HashMap<String, VCFRecord>) -> Vec<&'a VCFRecord> {
    let mut depth_sorted_records = Vec::new();

    for (_key, value) in id_to_record {
        depth_sorted_records.push(value);
    }

    depth_sorted_records.sort_by(|a, b| get_vcf_level(b).cmp(&get_vcf_level(a))); // reverse sort

    depth_sorted_records
}

// scan the levels (top-down) and write the output genotypes in id_to_genome.
// for level 0, the genotype is taken from pangenie (pos_to_pg_record)
// for anything else, it comes from the parent (and is looked up in the output table)
fn resolve_genotypes(id_to_record : &HashMap<String, VCFRecord>,
                     pos_to_pg_record : &HashMap<(String, u64), VCFRecord>,
                     levels : &Vec<&VCFRecord>,
                     id_to_genotype : &mut HashMap<String, Vec<u32>>) {
    
    for &vcf_record in levels {
        if get_vcf_level(vcf_record) == 0 {
            let chrom_bytes = vcf_record.chromosome.clone();
            let chrom_string = String::from_utf8(chrom_bytes).unwrap();
            println!("resolving {} {}", chrom_string, vcf_record.position);
            match pos_to_pg_record.get(&(chrom_string, vcf_record.position)) {
                Some(pg_vcf_record) => {
                    // dig out the first sample in the pangenie VCF
                    // todo: allow multiple samples (which should be pretty trivial)
                    let sample_bytes = &pg_vcf_record.header().samples()[0];
                    println!("finding sample {}", String::from_utf8(sample_bytes.to_vec()).unwrap());
                    match vcf_record.genotype(b"sample", b"GT") {
                        Some(gt) => {
                            let gt_bytes = gt[0].clone();
                            let gt_string = String::from_utf8(gt_bytes).unwrap();
                            println!("GIT STRING IS {}", gt_string);
                        },
                        None => panic!("couldnt find genotype"),
                    };
                    
                }
                None => panic!("Unable to find position {} in pangenie VCF", &vcf_record.position),
            };
        }
    }
}

fn main() -> Result<(), VCFError> {

	 let args: Vec<String> = env::args().collect();
	 let full_vcf_path = &args[1];
    let pg_vcf_path = &args[2];

    // index the full vcf from deconstruct (it must contain all the levels and annotations)
    println!("Indexing deconstruct VCF by ID: {}", full_vcf_path);
    let id_to_record = index_vcf_by_id(full_vcf_path).unwrap();
    // index the pangenie vcf.  its id's aren't consistent so we use coordinates instead
    println!("Indexing pangenie VCF by position: {}", pg_vcf_path);
    let pos_to_pg_record = index_vcf_by_pos(pg_vcf_path).unwrap();
    // in order to do a top-down traversal of the vcf, we sort it by decreasing level in the snarl
    // tree.  this information comes out of the LV info tag
    println!("Sorting by Level");
    let levels = sort_by_level(&id_to_record);
    
    // this is a map of resolved genotypes
    // todo: learn enough rust to be able to do it in place
    println!("Resolving genotypes");
    let mut id_to_genotype = HashMap::new();
    resolve_genotypes(&id_to_record, &pos_to_pg_record, &levels, &mut id_to_genotype);
    
    Ok(())
}
