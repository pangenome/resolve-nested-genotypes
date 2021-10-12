use std::collections::{HashMap, HashSet};

// represent snarl nesting relationships
// everything goes through integer ids to avoid refernces (which i don't seem to understand
// well enough to write event the simplest data structures)
pub struct SnarlForest {
    no_to_id : Vec<String>,
    id_to_no : HashMap<String, usize>,
    root_nos : HashSet<usize>,
    to_parent : HashMap<usize, usize>,
    to_children : HashMap<usize, Vec<usize>>,
}

impl SnarlForest {

    pub fn new() -> SnarlForest {
        SnarlForest {
            no_to_id : Vec::new(),
            id_to_no : HashMap::new(),
            root_nos : HashSet::new(),
            to_parent : HashMap::new(),
            to_children : HashMap::new(),
        }
    }
    
    // add node/parent to forest
    pub fn add_id(&mut self, id : &str, parent_id : Option<String>) -> usize {
        let no : usize;
        if self.id_to_no.contains_key(id) {
            no = *self.id_to_no.get(id).unwrap();
        } else {
            no = self.no_to_id.len();
            self.no_to_id.push(id.to_string());
            self.id_to_no.insert(id.to_string(), no);
        }
        if parent_id.is_some() {
            let par_no = self.add_id(&parent_id.unwrap(), None);
            if !self.to_parent.contains_key(&no) {
                self.to_parent.insert(no, par_no);
            }
            if !self.to_children.contains_key(&par_no) {
                self.to_children.insert(par_no, Vec::new());
            }
            self.to_children.get_mut(&par_no).unwrap().push(no);
        }
        no
    }

    pub fn calculate_roots(&mut self) {
        for no in 0..self.no_to_id.len() {
            if !self.to_parent.contains_key(&no) {
                self.root_nos.insert(no);
            }
        }
    }
}

