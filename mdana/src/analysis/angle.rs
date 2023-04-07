use super::super::struct_file::xtc::Coords;

pub fn angle(traj:&Vec<Coords>,idx1:usize,idx2:usize,idx3:usize)->Vec<f32>{
    let mut angvec=vec![];
    for i in 0..traj.len() {
        let a=((traj[i].coords[idx1][0]-traj[i].coords[idx2][0]).powi(2)+(traj[i].coords[idx1][1]-traj[i].coords[idx2][1]).powi(2)+(traj[i].coords[idx1][2]-traj[i].coords[idx2][2]).powi(2)).sqrt()/3.0;
        let b=((traj[i].coords[idx3][0]-traj[i].coords[idx2][0]).powi(2)+(traj[i].coords[idx3][1]-traj[i].coords[idx2][1]).powi(2)+(traj[i].coords[idx3][2]-traj[i].coords[idx2][2]).powi(2)).sqrt()/3.0;
        let c=((traj[i].coords[idx3][0]-traj[i].coords[idx1][0]).powi(2)+(traj[i].coords[idx3][1]-traj[i].coords[idx1][1]).powi(2)+(traj[i].coords[idx3][2]-traj[i].coords[idx1][2]).powi(2)).sqrt()/3.0;
        let cos=(a*a+b*b-c*c)/(2.0*a*b);
        let theta=cos.acos()*180.0/(std::f32::consts::PI);
        angvec.push(theta);
    }
    return angvec
}