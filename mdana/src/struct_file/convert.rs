pub struct Convert{
    pub data:Vec<String>,
}

impl Convert{
    pub fn vec_str_to_i32(&self)-> Vec<usize> {
        let mut result:Vec<usize>=vec![];
        for i in self.data.iter(){
           let input:usize=i.clone().trim().parse().expect("Please input a number vec");
           result.push(input);
        }
        return result
    }
}
