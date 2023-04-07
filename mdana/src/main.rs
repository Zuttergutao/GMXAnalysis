#![allow(dead_code, unused_imports, unused_mut, unused_variables)]
// #![warn(unused_assignments)]

mod analysis;
mod struct_file;

use crate::input::input;
use analysis::total::*;
use analysis::{angle, distance, input, mmpbsa, rmsd, xvg};
use std::io;
use std::time::{Duration, SystemTime};
use struct_file::gro_select::GroData;
use struct_file::xtc::Coords;
use struct_file::{convert, gro, gro_select, xtc};

use std::fs::File;
use std::io::{BufRead, BufReader};

use regex::Regex;
use std::collections::HashMap;
// todo
// map filter to calculate rmsd distance angle

fn main() {
    // let t1 = SystemTime::now();
    // let grodata= gro::gro("test/mdm.gro");
    // println!("{:?}",SystemTime::now().duration_since(t1).unwrap().as_secs());
    // let test=GroData{data:grodata}.select("name C1");
    // let filter=convert::Convert{data:test[0].clone()}.vec_str_to_i32();
    //     let m=xtc::loopxtc("test/mdm.xtc",0,100,1).unwrap();
    // let x=xtc::filterxtc("test/mdm.xtc",0,200,1,&filter).unwrap();
    //     println!("1-2 dist: {:?}",XTCdata{data:m}.distance(1, 2));
    // println!("1-2-3 angle: {:?}",angle::angle(&m,1,2,3));
    // println!("Rmsd: {:?}",rmsd::rmsd(&m,0));
    // println!("{:?}",select(&grodata, "resn LYS,ALA"));
    // println!("{:?}",select(&grodata, "resi 1,2,3:8"));
    // println!("{:?}",gro_select::select(&grodata, "resn PHE")[0]);
    // select(&grodata, "name C*");
    // xvg::readxvg("test/rmsf.xvg");
    // xvg::readxvg1("test/nhbdist.xvg",5);
    //    input();
    let file = File::open("md.out").unwrap();
    let lines = BufReader::new(file).lines();
    let mut molblockstart: Vec<usize> = vec![];
    // 盒子信息
    let mut bbox: Vec<Vec<f64>> = vec![vec![]];
    // 原子坐标
    let mut coord: Vec<Vec<f64>> = vec![vec![]];
    // 原子数
    let mut xlen: usize;
    // 体系block
    let mut blocks: HashMap<String, usize> = HashMap::new();
    // tpr文件中原子坐标开始行
    let mut xstart: usize;
    // tpr文件中原子坐标结束行
    let mut xend: usize;
    let mut tmp = 9999999;

    for (linenum, line) in lines.enumerate() {
        let mut m = line.unwrap();
        // 读取每个molblock起始行
        if m.starts_with("   molblock") {
            molblockstart.push(linenum);
        } else if m.starts_with("box (3x3):") {
            tmp = linenum;
        } else if m.starts_with("x (") {
            // 读取坐标
            let rexlen = Regex::new(r"\s\(([0-9]+)x").unwrap();
            // 读取原子数
            let xlen: usize = rexlen
                .captures(&m)
                .unwrap()
                .get(1)
                .map_or("", |m| m.as_str())
                .parse()
                .unwrap();
            let xstart = linenum + 2;
            let xend = xlen + xstart;
        } else if m.starts_with("   x[") {
            // 原子坐标
            let reatom = Regex::new(r"=\{([\s0-9a-zA-Ze\.\+\-,]+)\}").unwrap();
            let atom_vec: Vec<&str> = reatom
                .captures(&m)
                .unwrap()
                .get(1)
                .map_or("", |s| s.as_str())
                .trim()
                .split(",")
                .collect();
            let coord_vec: Vec<f64> = atom_vec.iter().map(|x| x.trim().parse().unwrap()).collect();
            coord.push(coord_vec);
        };
        // 读取盒子信息
        if (linenum == tmp + 1) | (linenum == tmp + 2) | (linenum == tmp + 3) {
            // 正则表达式提取盒子信息
            let re = Regex::new(r"=\{([,\s\+\-e\.0-9]+)\}").unwrap();
            // 提取盒子信息
            let re_vec: Vec<&str> = re
                .captures(&m)
                .unwrap()
                .get(1)
                .map_or("", |s| s.as_str())
                .trim()
                .split(",")
                .collect();
            // 处理Vec，将&str转化为f64
            let pro_vec: Vec<f64> = re_vec.iter().map(|x| x.trim().parse().unwrap()).collect();
            bbox.push(pro_vec);
        };
    }
    bbox.remove(0);
    coord.remove(0);
    println!("{:?}", coord[0]);
    println!("{:?}", bbox);
    println!("{:?}", molblockstart);
}
