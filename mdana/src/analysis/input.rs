use std::io;
use super::super::struct_file::{convert, gro, gro_select, xtc};
use super::super::XTCdata;

pub fn input() {
    let mut grovar = "空";
    let mut xtcvar = "空";
    let mut gropath: String;
    let mut xtcpath: String;
    'outer:loop {
        println!(
                "Please choose the function:\n-1: quit\n1: grofile: {}\n2: xtcfile: {}\n3: analysis",
        grovar, xtcvar
    );
        let mut input = String::new();
        io::stdin()
            .read_line(&mut input)
            .expect("Failed to read line");
        let input: i32 = input.trim().parse().expect("Please type a number!");
        match input {
            -1 => {
                println!("Enjoy the life!");
                break;
            }
            1 => loop {
                println!("Gro Function");
                println!("请输入gro文件路径，默认路径为当前目录");
                gropath = String::new();
                io::stdin()
                    .read_line(&mut gropath)
                    .expect("Failed to read path");
                if gropath.trim().get(1..2).unwrap() != ":" {
                    gropath = std::env::current_dir()
                        .unwrap()
                        .to_str()
                        .unwrap()
                        .to_owned()
                        + "\\"
                        + &gropath.trim().to_owned();
                }
                grovar = &gropath;
                if std::path::Path::new(grovar).exists() {
                    println!("Please choose the function:\n-1: quit\n1: grofile: {}\n2: xtcfile: {}\n3: analysis", grovar,xtcvar);
                    break;
                } else {
                    println!("Error!!! \"{}\" is not exist! \n", grovar);
                }
            },
            2 => loop {
                println!("XTC Function");
                println!("请输入xtc文件路径，默认路径为当前目录");
                xtcpath = String::new();
                io::stdin()
                    .read_line(&mut xtcpath)
                    .expect("Failed to read path");
                if xtcpath.trim().get(1..2).unwrap() != ":" {
                    xtcpath = std::env::current_dir()
                        .unwrap()
                        .to_str()
                        .unwrap()
                        .to_owned()
                        + "\\"
                        + &xtcpath.trim().to_owned();
                }
                xtcvar = &xtcpath;
                if std::path::Path::new(xtcvar).exists() {
                    println!("Please choose the function:\n-1: quit\n1: grofile: {}\n2: xtcfile: {}\n3: analysis", grovar,xtcvar);
                    break;
                } else {
                    println!("Error!!! \"{}\" is not exist! \n", xtcvar);
                }
            },
            3 => 'inner: loop {
                println!("Analysis Function");
                println!("-1 quit\n0 return \n1 Distance\n2 Angle\n3 RMSD\n4 Draw xvg");
                let mut input = String::new();
                io::stdin()
                    .read_line(&mut input)
                    .expect("Failed to read line");
                let input: i32 = input.trim().parse().expect("Please type a number!");
                match input {
                    -1 => break 'outer,
                    0 => {
                        println!("To do");
                        println!("-1 quit\n0 return \n1 Distance\n2 Angle\n3 RMSD\n4 Draw xvg");
                        continue 'inner;
                    }
                    1 => {
                        println!("select the atom 1");
                        let mut at1=String::new();
                        io::stdin().read_line(&mut at1).expect("Failed to read line");
                        let at1:usize=at1.trim().parse().expect("Please type a number");
                        println!("select the atom 2");
                        let mut at2=String::new();
                        io::stdin().read_line(&mut at2).expect("Failed to read line");
                        let at2:usize=at2.trim().parse().expect("Please type a number");
                        let m=xtc::loopxtc(xtcvar, 0, 100, 1).unwrap();
                        println!("{}-{} dist:{:?}",at1,at2,XTCdata{data:m}.distance(at1,at2));
                        continue 'inner;
                    }
                    2 => {
                        println!("select the atom 1");
                        let mut at1=String::new();
                        io::stdin().read_line(&mut at1).expect("Failed to read line");
                        let at1:usize=at1.trim().parse().expect("Please type a number");
                        println!("select the atom 2");
                        let mut at2=String::new();
                        io::stdin().read_line(&mut at2).expect("Failed to read line");
                        let at2:usize=at2.trim().parse().expect("Please type a number");
                        println!("select the atom 3");
                        let mut at3=String::new();
                        io::stdin().read_line(&mut at3).expect("Failed to read line");
                        let at3:usize=at3.trim().parse().expect("Please type a number");
                        let m=xtc::loopxtc(xtcvar, 0, 100, 1).unwrap();
                        println!("{}-{}-{} angle:{:?}",at1,at2,at3,XTCdata{data:m}.angle(at1, at2, at3));
                        continue 'inner;
                    }
                    3 => {
                        println!("select the reference frame");
                        let mut frame=String::new();
                        io::stdin().read_line(&mut frame).expect("Failed to read line");
                        let frame:usize=frame.trim().parse().expect("Please type a number");
                        let m=xtc::loopxtc(xtcvar, 0, 100, 1).unwrap();
                        println!("RMSD: {:?}",XTCdata{data:m}.rmsd(frame));
                        continue 'inner;
                    }
                    4 => {}
                    _ => {}
                }
            }

            _ => {
                continue;
            }
        }
    }
}