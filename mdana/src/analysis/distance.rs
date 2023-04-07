
use super::super::struct_file::xtc::Coords;


pub fn distance(traj: &Vec<Coords>,idx1:usize,idx2:usize)->Vec<f32>{
    let mut disvec=vec![];
    let mut x=traj[0].index.get(&idx1);
    // traj[0].index.get(&idx2);
    match x {
        Some(x)=> {
            let mut idx_1=*traj[0].index.get(&idx1).unwrap();
            let mut idx_2=*traj[0].index.get(&idx2).unwrap();
            for i in 0..traj.len() {
                let  tmp=((traj[i].coords[idx_1][0]-traj[i].coords[idx_2][0]).powi(2)+(traj[i].coords[idx_1][1]-traj[i].coords[idx_2][1]).powi(2)+(traj[i].coords[idx_1][2]-traj[i].coords[idx_2][2]).powi(2)).sqrt();
                disvec.push(tmp);
            }
            return disvec
        }
        None => {
            println!("*****************************");
            println!("Atom index is not in the traj, \nPlease check the atom index!");
            let mut vec =vec![];
            for key in traj[0].index.keys(){
                vec.push(key.clone());
            }
            vec.sort();
            println!("Atom index includes: {:?}",vec);
            println!("*****************************");
            return disvec
        }
    }

}