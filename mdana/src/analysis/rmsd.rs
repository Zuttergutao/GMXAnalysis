use super::super::struct_file::xtc::Coords;


pub fn rmsd(traj:&Vec<Coords>,rf:usize)->Vec<f32>{
    let mut rmsdvec=vec![];
    for i in 0..(traj.len()){
        let mut tmp:f32=0.0;
        for j in 0..(traj[rf].atoms as usize) {
            tmp += (traj[i].coords[j][0]-traj[rf].coords[j][0]).powi(2)+(traj[i].coords[j][1]-traj[rf].coords[j][1]).powi(2)+(traj[i].coords[j][2]-traj[rf].coords[j][2]).powi(2);
        }
        let rmsd=(tmp/(traj[rf].atoms as f32)).sqrt();
        rmsdvec.push(rmsd)
    }
    return rmsdvec
}