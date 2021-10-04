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
                    let id_bytes = vcf_record.id[0].clone();
                    let id_string = String::from_utf8(id_bytes).unwrap();
                    mymap.insert((id_string, vcf_record.position), vcf_record.clone());                    
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

fn main() -> Result<(), VCFError> {

	 let args: Vec<String> = env::args().collect();
	 let full_vcf_path = &args[1];
    let pg_vcf_path = &args[2];

    let id_to_record = index_vcf_by_id(full_vcf_path).unwrap();
    let pos_to_pg_recrd = index_vcf_by_pos(pg_vcf_path).unwrap();
    


    Ok(())
}
