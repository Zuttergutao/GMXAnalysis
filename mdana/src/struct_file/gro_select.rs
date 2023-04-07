/// To select atoms or residues for gro files
/// 
/// syntax for selection:
/// 
/// select residue index: `resi 1,2`  `resi 1,3,5:8`
/// 
/// select residue name: `resn LYS` `resn PHE`
/// 
/// select atom names: `name CA` `name CD*` `name CA,CB`
/// 
/// select atom index: `index 1` `index 1,2` `index 1:2`
/// return a Vec\<Vec\<String\>\> 
/// 
/// which contains \[\[atomidx\], \[resnidx\], \[resname\], \[atomname\], \[coords\], \[velocity\]\]
/// 
/// Warning!: 
/// 
/// No matter what the atomindex is in the gro file, the select function always starts at 0
/// ## Example
/// select C1 atoms in gro file
/// ```rust
/// let grodata= gro::gro("test/mdm.gro");
/// let test=GroData{data:grodata}.select("resi 1,2:8");
/// ```
pub  struct GroData {
   pub data:Vec<Vec<String>>,
}

impl GroData{
   pub fn select(&self,idx:&str)-> Vec<Vec<String>>{
      let sele_name:Vec<&str>=idx.split(" ").collect();
      let mut selet=vec![vec![]];
      let mut index=vec![];
      let mut resi=vec![];
      let mut resn=vec![];
      let mut atmn=vec![];
      let mut coord=vec![];
      let mut velocity=vec![];
      match sele_name[0] {
         "resi"=>{
            let mut tmp=vec![];
            for i in &sele_name[1].split(",").collect::<Vec<&str>>() {
               let atomidx=i.split(":").collect::<Vec<&str>>();
               if atomidx.len()==2{
                  let start:u32=atomidx[0].parse().unwrap() ;
                  let end:u32=atomidx[1].parse().unwrap();
                  if start > end {
                     panic!("Error input");
                  } else {
   
                     for i in start..end+1{
                        tmp.push(String::from(i.to_string()));
                     }
                  }
               } else {
                  tmp.push(String::from(i.to_string()));
               }
            }
            for (idx,res) in self.data[0].iter().enumerate(){
               for i in 0..(tmp.len()){
                  if tmp[i] == res.to_string() {
                     index.push(idx.to_string());
                     resi.push(self.data[0][idx].clone());
                     resn.push(self.data[1][idx].clone());
                     atmn.push(self.data[2][idx].clone());
                     coord.push(self.data[3][idx].clone());
                     velocity.push(self.data[4][idx].clone());
   
                  }
               }
            }
            selet.push(index);
            selet.push(resi);
            selet.push(resn);
            selet.push(atmn);
            selet.push(coord);
            selet.push(velocity);
            // remove the first null vec 
            // when push vec to  a vec, the first vec is empty
            selet.remove(0);
            return selet
         },
         "resn" => {
            let mut tmp=vec![];
            for i in &sele_name[1].split(",").collect::<Vec<&str>>() {
                  tmp.push(String::from(i.to_string()));
            }
            for (idx,res) in self.data[1].iter().enumerate(){
               for i in 0..(tmp.len()){
                  if tmp[i].to_uppercase() == res.to_string() {
                     index.push(idx.to_string());
                     resi.push(self.data[0][idx].clone());
                     resn.push(self.data[1][idx].clone());
                     atmn.push(self.data[2][idx].clone());
                     coord.push(self.data[3][idx].clone());
                     velocity.push(self.data[4][idx].clone());
   
                  }
               }
            }
            selet.push(index);
            selet.push(resi);
            selet.push(resn);
            selet.push(atmn);
            selet.push(coord);
            selet.push(velocity);
            // remove the first null vec 
            // when push vec to  a vec, the first vec is empty
            selet.remove(0);
            return selet
         },
         "name" => {
            let mut tmp=vec![];
            for i in &sele_name[1].split(",").collect::<Vec<&str>>() {
                  tmp.push(String::from(i.to_string()));
            }
            for (idx,res) in self.data[2].iter().enumerate(){
               for i in 0..(tmp.len()){
                  if tmp[i].ends_with("*") {
                     let idx_tmp=tmp[i].len()-1;
                     if res.len() > idx_tmp {
                        if &tmp[i][0..idx_tmp].to_uppercase() == &res[0..idx_tmp] {
                           index.push(idx.to_string());
                           resi.push(self.data[0][idx].clone());
                           resn.push(self.data[1][idx].clone());
                           atmn.push(self.data[2][idx].clone());
                           coord.push(self.data[3][idx].clone());
                           velocity.push(self.data[4][idx].clone());
                        }
                     }
                  }else{
                     if tmp[i].to_uppercase() == res.to_string() {
                        index.push(idx.to_string());
                        resi.push(self.data[0][idx].clone());
                        resn.push(self.data[1][idx].clone());
                        atmn.push(self.data[2][idx].clone());
                        coord.push(self.data[3][idx].clone());
                        velocity.push(self.data[4][idx].clone());
                     }
                  }
               }
            }
            selet.push(index);
            selet.push(resi);
            selet.push(resn);
            selet.push(atmn);
            selet.push(coord);
            selet.push(velocity);
            // remove the first null vec 
            // when push vec to  a vec, the first vec is empty
            selet.remove(0);
            return selet
         },
         "index"=>{
          let mut tmp=vec![];
          for i in &sele_name[1].split(",").collect::<Vec<&str>>() {
             let atomidx=i.split(":").collect::<Vec<&str>>();
             if atomidx.len()==2{
                let start:u32=atomidx[0].parse().unwrap() ;
                let end:u32=atomidx[1].parse().unwrap();
                if start > end {
                   panic!("Error input");
                } else {
  
                   for i in start..end+1{
                      tmp.push(String::from(i.to_string()));
                   }
                }
             } else {
                tmp.push(String::from(i.to_string()));
             }
          }
          for (idx,res) in self.data[0].iter().enumerate(){
             for i in 0..(tmp.len()){
                if tmp[i] == idx.to_string() {
                   index.push(idx.to_string());
                   resi.push(self.data[0][idx].clone());
                   resn.push(self.data[1][idx].clone());
                   atmn.push(self.data[2][idx].clone());
                   coord.push(self.data[3][idx].clone());
                   velocity.push(self.data[4][idx].clone());
  
                }
             }
          }
          selet.push(index);
          selet.push(resi);
          selet.push(resn);
          selet.push(atmn);
          selet.push(coord);
          selet.push(velocity);
          // remove the first null vec 
          // when push vec to  a vec, the first vec is empty
          selet.remove(0);
          return selet
          },
         _=>{
            return selet
         },
      }
   }


}


// pub fn select(gro:&Vec<Vec<String>>,idx:&str)-> Vec<Vec<String>>{
//    let sele_name:Vec<&str>=idx.split(" ").collect();
//    let mut selet=vec![vec![]];
//    let mut index=vec![];
//    let mut resi=vec![];
//    let mut resn=vec![];
//    let mut atmn=vec![];
//    let mut coord=vec![];
//    let mut velocity=vec![];
//    match sele_name[0] {
//       "resi"=>{
//          let mut tmp=vec![];
//          for i in &sele_name[1].split(",").collect::<Vec<&str>>() {
//             let atomidx=i.split(":").collect::<Vec<&str>>();
//             if atomidx.len()==2{
//                let start:u32=atomidx[0].parse().unwrap() ;
//                let end:u32=atomidx[1].parse().unwrap();
//                if start > end {
//                   panic!("Error input");
//                } else {

//                   for i in start..end+1{
//                      tmp.push(String::from(i.to_string()));
//                   }
//                }
//             } else {
//                tmp.push(String::from(i.to_string()));
//             }
//          }
//          for (idx,res) in gro[0].iter().enumerate(){
//             for i in 0..(tmp.len()){
//                if tmp[i] == res.to_string() {
//                   index.push(idx.to_string());
//                   resi.push(gro[0][idx].clone());
//                   resn.push(gro[1][idx].clone());
//                   atmn.push(gro[2][idx].clone());
//                   coord.push(gro[3][idx].clone());
//                   velocity.push(gro[4][idx].clone());

//                }
//             }
//          }
//          selet.push(index);
//          selet.push(resi);
//          selet.push(resn);
//          selet.push(atmn);
//          selet.push(coord);
//          selet.push(velocity);
//          // remove the first null vec 
//          // when push vec to  a vec, the first vec is empty
//          selet.remove(0);
//          return selet
//       },
//       "resn" => {
//          let mut tmp=vec![];
//          for i in &sele_name[1].split(",").collect::<Vec<&str>>() {
//                tmp.push(String::from(i.to_string()));
//          }
//          for (idx,res) in gro[1].iter().enumerate(){
//             for i in 0..(tmp.len()){
//                if tmp[i].to_uppercase() == res.to_string() {
//                   index.push(idx.to_string());
//                   resi.push(gro[0][idx].clone());
//                   resn.push(gro[1][idx].clone());
//                   atmn.push(gro[2][idx].clone());
//                   coord.push(gro[3][idx].clone());
//                   velocity.push(gro[4][idx].clone());

//                }
//             }
//          }
//          selet.push(index);
//          selet.push(resi);
//          selet.push(resn);
//          selet.push(atmn);
//          selet.push(coord);
//          selet.push(velocity);
//          // remove the first null vec 
//          // when push vec to  a vec, the first vec is empty
//          selet.remove(0);
//          return selet
//       },
//       "name" => {
//          let mut tmp=vec![];
//          for i in &sele_name[1].split(",").collect::<Vec<&str>>() {
//                tmp.push(String::from(i.to_string()));
//          }
//          for (idx,res) in gro[2].iter().enumerate(){
//             for i in 0..(tmp.len()){
//                if tmp[i].ends_with("*") {
//                   let idx_tmp=tmp[i].len()-1;
//                   if res.len() > idx_tmp {
//                      if &tmp[i][0..idx_tmp].to_uppercase() == &res[0..idx_tmp] {
//                         index.push(idx.to_string());
//                         resi.push(gro[0][idx].clone());
//                         resn.push(gro[1][idx].clone());
//                         atmn.push(gro[2][idx].clone());
//                         coord.push(gro[3][idx].clone());
//                         velocity.push(gro[4][idx].clone());
//                      }
//                   }
//                }else{
//                   if tmp[i].to_uppercase() == res.to_string() {
//                      index.push(idx.to_string());
//                      resi.push(gro[0][idx].clone());
//                      resn.push(gro[1][idx].clone());
//                      atmn.push(gro[2][idx].clone());
//                      coord.push(gro[3][idx].clone());
//                      velocity.push(gro[4][idx].clone());
//                   }
//                }
//             }
//          }
//          selet.push(index);
//          selet.push(resi);
//          selet.push(resn);
//          selet.push(atmn);
//          selet.push(coord);
//          selet.push(velocity);
//          // remove the first null vec 
//          // when push vec to  a vec, the first vec is empty
//          selet.remove(0);
//          return selet
//       },
//       "index"=>{
//        let mut tmp=vec![];
//        for i in &sele_name[1].split(",").collect::<Vec<&str>>() {
//           let atomidx=i.split(":").collect::<Vec<&str>>();
//           if atomidx.len()==2{
//              let start:u32=atomidx[0].parse().unwrap() ;
//              let end:u32=atomidx[1].parse().unwrap();
//              if start > end {
//                 panic!("Error input");
//              } else {

//                 for i in start..end+1{
//                    tmp.push(String::from(i.to_string()));
//                 }
//              }
//           } else {
//              tmp.push(String::from(i.to_string()));
//           }
//        }
//        for (idx,res) in gro[0].iter().enumerate(){
//           for i in 0..(tmp.len()){
//              if tmp[i] == idx.to_string() {
//                 index.push(idx.to_string());
//                 resi.push(gro[0][idx].clone());
//                 resn.push(gro[1][idx].clone());
//                 atmn.push(gro[2][idx].clone());
//                 coord.push(gro[3][idx].clone());
//                 velocity.push(gro[4][idx].clone());

//              }
//           }
//        }
//        selet.push(index);
//        selet.push(resi);
//        selet.push(resn);
//        selet.push(atmn);
//        selet.push(coord);
//        selet.push(velocity);
//        // remove the first null vec 
//        // when push vec to  a vec, the first vec is empty
//        selet.remove(0);
//        return selet
//        },
//       _=>{
//          return selet
//       },
//    }

// }

 