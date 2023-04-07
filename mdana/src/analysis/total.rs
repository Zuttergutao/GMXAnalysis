use super::super::struct_file::xtc::Coords;

pub struct XTCdata{
    pub data:Vec<Coords>,
}


impl XTCdata {
    pub fn distance(&self,idx1:usize,idx2:usize)->Vec<f32>{
        let mut disvec=vec![];
        for i in 0..self.data.len() {
            let  tmp=((self.data[i].coords[idx1][0]-self.data[i].coords[idx2][0]).powi(2)+(self.data[i].coords[idx1][1]-self.data[i].coords[idx2][1]).powi(2)+(self.data[i].coords[idx1][2]-self.data[i].coords[idx2][2]).powi(2)).sqrt();
            disvec.push(tmp);
        }
        let mut avgdis:f32=0.0;
        for i in disvec.iter() {
            avgdis = avgdis+i;
        }
        avgdis=avgdis / disvec.len() as f32;
        println!("The average distance between {} and {} is {} nm",idx1,idx2,avgdis);
        return disvec
    }

    pub fn angle(&self,idx1:usize,idx2:usize,idx3:usize)->Vec<f32>{
        let mut angvec=vec![];
        for i in 0..self.data.len() {
            let a=((self.data[i].coords[idx1][0]-self.data[i].coords[idx2][0]).powi(2)+(self.data[i].coords[idx1][1]-self.data[i].coords[idx2][1]).powi(2)+(self.data[i].coords[idx1][2]-self.data[i].coords[idx2][2]).powi(2)).sqrt()/3.0;
            let b=((self.data[i].coords[idx3][0]-self.data[i].coords[idx2][0]).powi(2)+(self.data[i].coords[idx3][1]-self.data[i].coords[idx2][1]).powi(2)+(self.data[i].coords[idx3][2]-self.data[i].coords[idx2][2]).powi(2)).sqrt()/3.0;
            let c=((self.data[i].coords[idx3][0]-self.data[i].coords[idx1][0]).powi(2)+(self.data[i].coords[idx3][1]-self.data[i].coords[idx1][1]).powi(2)+(self.data[i].coords[idx3][2]-self.data[i].coords[idx1][2]).powi(2)).sqrt()/3.0;
            let cos=(a*a+b*b-c*c)/(2.0*a*b);
            let theta=cos.acos()*180.0/(std::f32::consts::PI);
            angvec.push(theta);
        }
        let mut avgang:f32=0.0;
        for i in angvec.iter() {
            avgang = avgang+i;
        }
        avgang=avgang / angvec.len() as f32;
        println!("The average angle between {}, {} and {} is {} degrees",idx1,idx2,idx3,avgang);
        return angvec
    }

    pub fn rmsd(&self,rf:usize)->Vec<f32>{
        let mut rmsdvec=vec![];
        for i in 0..(self.data.len()){
            let mut tmp:f32=0.0;
            for j in 0..(self.data[rf].atoms as usize) {
                tmp += (self.data[i].coords[j][0]-self.data[rf].coords[j][0]).powi(2)+(self.data[i].coords[j][1]-self.data[rf].coords[j][1]).powi(2)+(self.data[i].coords[j][2]-self.data[rf].coords[j][2]).powi(2);
            }
            let rmsd=(tmp/(self.data[rf].atoms as f32)).sqrt();
            rmsdvec.push(rmsd)
        }
        return rmsdvec
    }
}