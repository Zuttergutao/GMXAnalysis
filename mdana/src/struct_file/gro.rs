
use std::fs::File;
use std::io::{BufRead, BufReader};


/// Read Gro file 
/// 
/// \[0\] resindex 
/// 
/// \[1\] resname 
/// 
/// \[2\] atomname 
/// 
/// \[3\] coord 
/// 
/// \[4\] velocity 
/// 
/// ## Example
/// 
/// grodata=gro(&gro_path) 
pub fn gro(gro_path:&str)-> Vec<Vec<String>>{
    let mut resn=vec![];
    let mut resi=vec![];
    let mut atmn=vec![];
    let mut coord=vec![];
    let mut velocity=vec![];
    let mut total=vec![];
    let file=File::open(gro_path).unwrap();
    let lines=BufReader::new(file).lines();

    for (linenum,line) in lines.enumerate() {
        match linenum {
            0 => continue,
            1 => continue,
            _ => {
                if let Ok(data) =line {
                    if data.len() > 30 {
                        resi.push(data[0..5].trim().to_string());
                        resn.push(data[5..10].trim().to_string());
                        atmn.push(data[10..15].trim().to_string());
                        coord.push(data[20..44].trim().to_string());
                        velocity.push(data[44..68].trim().to_string())
                    } else {
                        continue
                    }
                }
            }
        }
    }
    total.push(resi);
    total.push(resn);
    total.push(atmn);
    total.push(coord);
    total.push(velocity);
    return total
}



     


// #[derive(Debug)]
// pub struct GroGeom(pub u64,pub String, pub String,pub f64,pub f64,pub f64);


// pub fn gro(gro_path:&str) -> HashMap<String,GroGeom> {
//     let file=File::open(gro_path).unwrap();
//     let lines=BufReader::new(file).lines();
//     let mut hashmap: HashMap<String, GroGeom>=HashMap::new();

//     for (linenum,line) in lines.enumerate() {
//         match linenum {
//             0 => continue,
//             1 => continue,
//             _ => {
//                 if let Ok(data) =line {
//                     if data.len() > 30 {
//                         let gro_data= GroGeom(data[0..5].trim().parse().unwrap(),data[5..10].trim().to_string(),data[10..15].trim().to_string(),data[20..28].trim().parse().unwrap(),data[28..36].trim().parse().unwrap(),data[36..44].trim().parse().unwrap());
//                         hashmap.insert(data[15..20].trim().parse().unwrap(), gro_data);
//                     } else {
//                         continue
//                     }
//                 }
//             }
//         }
//     }
//     return hashmap
// }

    
    



    




    
    



    
