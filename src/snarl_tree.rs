use std::collections::{HashMap, HashSet, VecDeque};

// represent snarl nesting relationships
// everything goes through integer ids to avoid refernces (which i don't seem to understand
// well enough to write event the simplest data structures)

// todo: no_to_id not needed
pub struct SnarlForest {
    no_to_id : Vec<String>,
    id_to_no : HashMap<String, usize>,
    root_nos : HashSet<usize>,
    to_parent : HashMap<usize, usize>,
    to_children : HashMap<usize, Vec<usize>>,
    no_to_hprc_id : HashMap<usize, String>,
}

impl SnarlForest {

    pub fn new() -> SnarlForest {
        SnarlForest {
            no_to_id : Vec::new(),
            id_to_no : HashMap::new(),
            root_nos : HashSet::new(),
            to_parent : HashMap::new(),
            to_children : HashMap::new(),
            no_to_hprc_id : HashMap::new(),
        }
    }
    
    // add node/parent to forest
    pub fn add_id(&mut self, id : &str, hprc_id : Option<String>, parent_id : Option<String>) -> usize {
        let no : usize;
        if self.id_to_no.contains_key(id) {
            no = *self.id_to_no.get(id).unwrap();        
        } else {
            no = self.no_to_id.len();
            self.no_to_id.push(id.to_string());
            self.id_to_no.insert(id.to_string(), no);
            if hprc_id.is_some() {
                self.no_to_hprc_id.insert(no, hprc_id.unwrap());
            }
        }
        if parent_id.is_some() {
            let par_no = self.add_id(&parent_id.unwrap(), None, None);
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

    // run once after adding everything to compute all the roots
    pub fn calculate_roots(&mut self) {
        for no in 0..self.no_to_id.len() {
            if !self.to_parent.contains_key(&no) {
                self.root_nos.insert(no);
            }
        }
    }
    
    // numeric to string id
    pub fn get_id(&self, no : usize) -> &str {
        &self.no_to_id[no]
    }

    // string to numeric id
    pub fn get_no(&self, id : &str) -> usize {
        *self.id_to_no.get(id).unwrap()
    }

    // return all leaves below and including node
    pub fn leaves_below(&self, id : &str) -> Vec<usize> {
        let mut leaves : Vec<usize> = Vec::new();

        let mut queue : VecDeque<usize> = VecDeque::new();
        queue.push_back(self.get_no(id));

        while queue.len() > 0 {
            let no = queue.pop_front().unwrap();
            match self.to_children.get(&self.get_no(id)) {
                Some(children) => {
                    for child in children {
                        queue.push_back(*child);
                    }
                },
                None => {
                    leaves.push(no);
                }
            };
        
        }
        leaves                
    }

    // leaf : id
    // not a leaf: colon-separated list of leaves below
    pub fn get_hprc_id(&self, id : &str) -> String {
        let no = self.id_to_no.get(id).expect("Could not find ID");
        let leaf_nos = self.leaves_below(id);
        let mut hprc_id = String::new();
        for i in 0..leaf_nos.len() {
            hprc_id.push_str(self.no_to_hprc_id.get(no).unwrap());
            if i < leaf_nos.len() - 1 {
                hprc_id.push(':');
            }
        }
        hprc_id
    }
}

