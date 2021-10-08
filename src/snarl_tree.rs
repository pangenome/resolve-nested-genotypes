use std::collections::{HashMap};
use std::rc::Rc;
use std::cell::{RefCell, RefMut};
use std::borrow::BorrowMut;

// tree node contains just a snarl id
struct SnarlTreeNode<'a> {
    id : &'a str,
    parent : Option<&'a RefCell<SnarlTreeNode<'a>>>,
    children : Vec<&'a RefCell<SnarlTreeNode<'a>>>,
}

// forest nodes are owned by the hashma which indexes them by name
// we wrap in RefCell to force rust to give us mutable access to references
// otherwise it's impossible (in my hands at least) to have multiple references
// and still want to to do updates. 
pub struct SnarlForest<'a> {
    roots : Vec<&'a RefCell<SnarlTreeNode<'a>>>,
    idx : HashMap<&'a str, RefCell<SnarlTreeNode<'a>>>,
}

impl<'a> SnarlForest<'a> {

    pub fn has_node(&'a self, id : &'a str) -> bool {
        return self.idx.contains_key(id);
    }

    pub fn get_node(&'a self, id : &'a str) -> &'a RefCell<SnarlTreeNode<'a>> {
        return self.idx.get(id).unwrap();
    }

    pub fn get_parent_id(&'a self, id : &'a str) -> Option<&'a str> {
        assert_eq!(self.has_node(id), true);
        match self.get_node(id).borrow().parent {
            Some(parent) => Some(parent.borrow().id),
            None => None,
        }
    }

    // todo: would be nice to provide iterator but am a bit too rusted out atm
    pub fn get_child_ids(&'a self, id : &'a str) -> Vec<&'a str> {
        assert_eq!(self.has_node(id), true);
        let mut children : Vec<&'a str> = Vec::new();
        for child in &self.get_node(id).borrow().children {
            children.push(child.borrow().id);            
        }
        children
    }

    // todo: ditto above
    pub fn get_root_ids(&'a self) -> Vec<&'a str> {
        let mut roots : Vec<&'a str> = Vec::new();
        for root in &self.roots {
            roots.push(root.borrow().id);
        }
        roots
    }
    
    // make a new node and add it to the forest (but not attached to anything)
    pub fn add_node(&'a mut self, id : &'a str) {
        let new_node = RefCell::new(SnarlTreeNode {
            id : id,
            parent : None,
            children : vec![]
        });
        self.idx.insert(id, new_node);
    }

    // set a parent
    pub fn set_parent(&'a self, child : &'a RefCell<SnarlTreeNode<'a>>, parent : &'a RefCell<SnarlTreeNode<'a>>) {
        child.borrow_mut().parent = Some(parent);
        parent.borrow_mut().children.push(child);
    }

    // figure out roots.  run once, after all nodes added and set parented.
    pub fn find_roots(&'a mut self) {
        for node in self.idx.values() {
            if node.borrow().parent.is_none() {
                self.roots.push(node);
            }
        }
    }
 
}
