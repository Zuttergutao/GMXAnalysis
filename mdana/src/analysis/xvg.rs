use plotly::common::Mode;
use plotly::layout::Legend;
use plotly::{ImageFormat, Plot, Scatter};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};

#[derive(Debug)]
struct XvgContents {
    title: String,
    xlabel: String,
    ylabel: String,
    xdata: Vec<f64>,
    ydata: Vec<f64>,
}

pub fn readxvg(xvg_path: &str) {
    let mut title: String = String::new();
    let mut xaxis: String = String::new();
    let mut yaxis: String = String::new();
    let mut xdata = Vec::new();
    let mut ydata = Vec::new();

    let mut file = File::open(xvg_path).unwrap();
    let mut lines = BufReader::new(file).lines();
    for (linenum, mut line) in lines.enumerate() {
        let mut m = line.unwrap();
        let mut s = &m.trim().get(0..1).unwrap().to_string();
        match &s as &str {
            "#" => continue,
            "@" => {
                let mut atinfo = &m
                    .trim_start_matches("@")
                    .trim()
                    .get(0..5)
                    .unwrap()
                    .to_string();
                match &atinfo as &str {
                    "title" => {
                        title = m
                            .trim_start_matches("@")
                            .trim()
                            .trim_start_matches("title")
                            .trim()
                            .to_string();
                    }
                    "xaxis" => {
                        xaxis = m
                            .trim_start_matches("@")
                            .trim()
                            .trim_start_matches("xaxis")
                            .trim()
                            .to_string();
                    }
                    "yaxis" => {
                        yaxis = m
                            .trim_start_matches("@")
                            .trim()
                            .trim_start_matches("yaxis")
                            .trim()
                            .to_string();
                    }
                    _ => continue,
                }
            }
            _ => {
                //  let data=line;
                for (idx, i) in m.split_ascii_whitespace().into_iter().enumerate() {
                    match idx {
                        0 => xdata.push(i.to_owned().parse().unwrap()),
                        1 => ydata.push(i.to_owned().parse().unwrap()),
                        _ => continue,
                    }
                }
            }
        }
    }
    let xvg = XvgContents {
        title: title,
        xlabel: xaxis,
        ylabel: yaxis,
        xdata: xdata,
        ydata: ydata,
    };

    let mut plot = Plot::new();
    let trace = Scatter::new(xvg.xdata, xvg.ydata);
    plot.add_trace(trace);

    // plot.write_image("out.png", ImageFormat::PNG, 800, 600, 1.0);
    // 输出到html
    plot.write_html("out.html");
}

#[derive(Debug)]
struct XvgContents1 {
    title: String,
    xlabel: String,
    ylabel: String,
    xdata: Vec<f64>,
    ydata: Vec<Vec<f64>>,
    legend: Vec<String>,
}

pub fn readxvg1(xvg_path: &str, dataset: usize) {
    let mut title: String = String::new();
    let mut xaxis: String = String::new();
    let mut yaxis: String = String::new();
    let mut legend = Vec::new();
    let mut legendname = Vec::new();
    let mut xdata = Vec::new();
    let mut ydata: Vec<Vec<f64>> = Vec::new();
    let mut data: Vec<Vec<f64>> = vec![Vec::new(); dataset];
    let mut file = File::open(xvg_path).unwrap();
    let mut lines = BufReader::new(file).lines();
    for (linenum, mut line) in lines.enumerate() {
        let mut m = line.unwrap();
        let mut s = &m.trim().get(0..1).unwrap().to_string();
        match &s as &str {
            "#" => continue,
            "@" => {
                let mut atinfo = &m
                    .trim_start_matches("@")
                    .trim()
                    .get(0..5)
                    .unwrap()
                    .to_string();
                match &atinfo as &str {
                    "title" => {
                        title = m
                            .trim_start_matches("@")
                            .trim()
                            .trim_start_matches("title")
                            .trim()
                            .to_string();
                    }
                    "xaxis" => {
                        xaxis = m
                            .trim_start_matches("@")
                            .trim()
                            .trim_start_matches("xaxis")
                            .trim()
                            .to_string();
                    }
                    "yaxis" => {
                        yaxis = m
                            .trim_start_matches("@")
                            .trim()
                            .trim_start_matches("yaxis")
                            .trim()
                            .to_string();
                    }
                    _ => {
                        if atinfo.contains("s") {
                            let mut tmp = m
                                .trim_start_matches("@")
                                .trim()
                                .trim_start_matches("yaxis")
                                .trim()
                                .get(0..3)
                                .unwrap()
                                .trim();
                            legend.push(tmp.to_owned());
                            legendname.push(
                                m.trim_start_matches("@")
                                    .trim()
                                    .trim_start_matches(tmp)
                                    .trim()
                                    .trim_start_matches("legend")
                                    .trim()
                                    .to_owned(),
                            );
                        }
                    }
                }
            }
            _ => {
                //  let data=line;
                for (idx, i) in m.split_ascii_whitespace().into_iter().enumerate() {
                    if idx == 0 {
                        xdata.push(i.to_owned().parse().unwrap());
                    } else {
                        data[idx - 1].push(i.to_owned().parse().unwrap());
                    }
                }
            }
        }
    }
    let xvg = XvgContents1 {
        title: title,
        xlabel: xaxis,
        ylabel: yaxis,
        xdata: xdata,
        ydata: data,
        legend: legendname,
    };
    println!("{:?}", xvg.legend)

    // let mut plot = Plot::new();
    // let trace = Scatter::new(xvg.xdata,xvg.ydata);
    // plot.add_trace(trace);

    // // plot.write_image("out.png", ImageFormat::PNG, 800, 600, 1.0);
    // // 输出到html
    // plot.write_html("out.html");
}
