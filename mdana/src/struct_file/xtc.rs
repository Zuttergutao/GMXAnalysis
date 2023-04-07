// use polars::export::arrow::compute::filter;
use std::collections::HashMap;
use std::time::{Duration, SystemTime};
use xdrfile::*;

#[derive(Debug)]
pub struct Coords {
    pub frame: f32,
    pub atoms: f32,
    pub index: HashMap<usize, usize>,
    pub time: f32,
    pub step: f32,
    pub coords: Vec<[f32; 3]>,
    pub box_vec: [[f32; 3]; 3],
}

pub struct Coords1 {
    pub frame: f32,
    pub atoms: f32,
    pub time: f32,
    pub step: f32,
    pub coords: Vec<[f32; 3]>,
    pub box_vec: [[f32; 3]; 3],
}



/// Read xtc file looply
///
/// Steps, begin frame, end frame can be set
///
/// Return a Vec contains frame, atoms, time, step, coords, box_vec
///
/// ## Example:
/// Read the second atomic info  in a step of 10 in the xtc file
///
/// atom index starts from 0
/// ```rust
/// let filter = vec![1];
/// let test = xtc::loopxtc("test/md.xtc",0,100,10);
/// ```
///
pub fn loopxtc(
    xtcfile: &str,
    begin: usize,
    end: usize,
    step: usize,
) -> Result<Vec<Coords>, Box<dyn std::error::Error>> {
    // 循环读取轨迹
    let mut xtc_trj = XTCTrajectory::open_read(xtcfile)?;
    let mut frame = Frame::with_len(xtc_trj.get_num_atoms()?);
    let mut vec = vec![];
    let mut count = 0;
    let mut map: HashMap<_, _> = HashMap::new();
    for v in xtc_trj.get_num_atoms().iter() {
        map.insert(*v, *v);
    }
    while let Ok(()) = xtc_trj.read(&mut frame) {
        if count >= begin && count <= end {
            match count % step {
                0 => {
                    let data = Coords {
                        frame: count as f32,
                        atoms: frame.len() as f32,
                        time: frame.time / 100.0,
                        step: frame.step as f32,
                        index: map.clone(),
                        coords: frame.coords.clone(),
                        box_vec: frame.box_vector,
                    };
                    vec.push(data);
                    count = count + 1;
                }
                _ => {
                    count = count + 1;
                    continue;
                }
            }
        } else {
            count = count + 1;
        }
    }
    let tframes = vec.len();
    let tbegin = vec[0].time;
    let tend = (vec[1].time - vec[0].time) * ((vec.len() - 1) as f32);
    let dt = vec[1].time - vec[0].time;
    println!(
        "Total frames read is {} from {} ps to {} ps  with dt {} ps",
        tframes, tbegin, tend, dt
    );
    Ok(vec)
}



/// Read filtered xtc file looply
///
/// Steps, begin frame, end frame, atoms filter can be set
///
/// Return a Vec contains frame, atoms, time, step, coords, box_vec
///
/// ## Example:
/// Read the second atomic info  in a step of 10 in the xtc file
///
/// atom index starts from 0
/// ```rust
/// let filter = vec![1];
/// let test = xtc::filterxtc("test/md.xtc",10,&filter);
/// ```
///
pub fn filterxtc(
    xtcfile: &str,
    begin: usize,
    end: usize,
    step: usize,
    filter: &Vec<usize>,
) -> Result<Vec<Coords>, Box<dyn std::error::Error>> {
    // 循环读取轨迹
    let mut xtc_trj = XTCTrajectory::open_read(xtcfile)?;
    let mut frame = Frame::with_len(xtc_trj.get_num_atoms()?);
    let mut vec = vec![];
    let mut count = 0;
    let mut map: HashMap<_, _> = HashMap::new();
    for (k, v) in filter.iter().enumerate() {
        map.insert(*v, k);
    }
    while let Ok(()) = xtc_trj.read(&mut frame) {
        let mut new_frame=&mut frame.clone();
        new_frame.filter_coords(filter);
        if count >= begin && count <= end {
            match count % step {
                0 => {
                    let data = Coords {
                        frame: count as f32,
                        atoms: new_frame.len() as f32,
                        index: map.clone(),
                        time: new_frame.time / 100.0,
                        step: new_frame.step as f32,
                        coords: new_frame.coords.clone(),
                        box_vec: new_frame.box_vector,
                    };
                    // println!("count: {:?}",count);
                    // println!("len: {:?}",new_frame.len());
                    // println!("time: {:?}",frame.time);
                    // println!("coords: {:?}",new_frame.coords.clone());
                    // println!("box: {:?}",new_frame.box_vector);
                    vec.push(data);
                    count = count + 1;
                }
                _ => {
                    count = count + 1;
                    continue;
                }
            }
        } else {
            count = count + 1;
        }
    }
    let tframes = vec.len();
    let tbegin = vec[0].time;
    let tend = (vec[1].time - vec[0].time) * ((vec.len() - 1) as f32);
    let dt = vec[1].time - vec[0].time;
    println!(
        "Total frames read is {} from {} ps to {} ps  with dt {} ps",
        tframes, tbegin, tend, dt
    );
    Ok(vec)
}




// pub fn loopxtc1(xtcfile:&str,begin:usize,end:usize,step:usize)-> Result<Vec<Coords1>> {
//     // 循环读取轨迹
//     let  trj = XTCTrajectory::open_read(xtcfile)?;
//     let mut vec=vec![];
//     for (idx,result) in trj.into_iter().enumerate(){
//         if idx >= begin && idx <= end{
//             match idx%step {
//                 0 => {
//                     let frame=result?;
//                     let data=Coords1 {
//                         frame:idx as f32,
//                         atoms:frame.len() as f32,
//                         time:frame.time/100.0,
//                         step:frame.step as f32,
//                         coords:frame.coords.clone(),
//                         box_vec:frame.box_vector,
//                     };
//                     vec.push(data)
//                 },
//                 _ => continue
//             }
//         }
//     };
//     let tframes=vec.len();
//     let tbegin=vec[0].time;
//     let tend=(vec[1].time-vec[0].time)*((vec.len()-1) as f32 );
//     let dt=vec[1].time-vec[0].time;
//     println!("Total frames read is {} from {} ps to {} ps  with dt {} ps", tframes,tbegin,tend,dt );
//     return Ok(vec)
// }

// 读取第一帧轨迹
// pub fn xtc(xtcfile:&str) -> Result<Coords> {
//     // 读取轨迹
//     let mut trj = XTCTrajectory::open_read(xtcfile)?;
//     let num_atoms = trj.get_num_atoms()?;
//     let mut frame = Frame::with_len(num_atoms);
//     trj.read(&mut frame)?;

//     let data=Coords {
//         atoms:frame.len() as f32,
//         time:frame.time as f32,
//         step:frame.step as f32,
//         // coords:frame.coords,
//         box_vec:frame.box_vector,
//     };
//     return Ok(data)
// }

// pub fn loopxtc(xtcfile:&str,begin:usize,end:usize,step:usize,filter:&Vec<usize>)->Result<Vec<Coords>, Box<dyn std::error::Error>>{
//     // 循环读取轨迹
//     let mut xtc_trj = XTCTrajectory::open_read(xtcfile)?;
//     let mut frame=Frame::with_len(xtc_trj.get_num_atoms()?);
//     let mut vec=vec![];
//     let mut count=0;
//     let mut map:HashMap<_, _>=HashMap::new();
//     for (k,v) in filter.iter().enumerate(){
//         map.insert(*v, k);
//     }
//     while let Ok(()) = xtc_trj.read(&mut frame){
//         // This step consume a lot of time
//         println!("{:?}",&mut frame.filter_coords(&[1,2]));
//         if count >= begin && count <= end {
//             match count%step {
//                 0=> {
//                     let data=Coords {
//                         frame:count as f32,
//                         atoms:frame.len() as f32,
//                         index:map.clone() ,
//                         time:frame.time/100.0,
//                         step:frame.step as f32,
//                         coords:frame.coords.clone(),
//                         box_vec:frame.box_vector,
//                     };
//                     // println!("count: {:?}",count);
//                     // println!("len: {:?}",new_frame.len());
//                     // println!("time: {:?}",new_frame.time);
//                     // println!("coords: {:?}",new_frame.coords.clone());
//                     // println!("box: {:?}",new_frame.box_vector);
//                     vec.push(data);
//                     count=count+1;
//                 },
//                 _ => {
//                     count=count+1;
//                     continue
//                 }
//             }
//         } else{
//             count=count+1;
//         }
//     }
//     let tframes=vec.len();
//     let tbegin=vec[0].time;
//     let tend=(vec[1].time-vec[0].time)*((vec.len()-1) as f32 );
//     let dt=vec[1].time-vec[0].time;
//     println!("Total frames read is {} from {} ps to {} ps  with dt {} ps", tframes,tbegin,tend,dt );
//     Ok(vec)
// }
