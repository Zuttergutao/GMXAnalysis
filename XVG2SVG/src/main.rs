use quill::*;
use quill::Plot as QuillPlot;
use std::fs::File;
use std::io::{BufRead, BufReader};
use clap::{Arg,App};
use textplots::{Chart, Shape,Plot};
use terminal_size::{Width, Height, terminal_size};

// https://github.com/Ryan-D-Gast/quill
#[derive(Debug)]
struct XvgMetadata {
    title: Option<String>,
    x_label: Option<String>,
    y_label: Option<String>,
}

#[derive(Debug)]
struct XvgData {
    metadata: XvgMetadata,
    data: Vec<Vec<f64>>,
}

fn parse_xvg_metadata_line(line: &str) -> Option<(String, String)> {

    if let Some(pos) = line.find('"') {
        let key_part = line[1..pos].trim(); 
        let rest = &line[pos+1..];
        if let Some(end_quote_pos) = rest.find('"') {
            let value = &rest[..end_quote_pos];
            return Some((key_part.to_string(), value.to_string()));
        }
    }
    None
}

fn read_xvg(path: &str) -> std::io::Result<XvgData> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut metadata = XvgMetadata {
        title: None,
        x_label: None,
        y_label: None,
    };
    let mut data_columns: Vec<Vec<f64>> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let line_trimmed = line.trim();

        if line_trimmed.is_empty() || line_trimmed.starts_with('#') {
            continue;
        } else if line_trimmed.starts_with('@') {
            if let Some((key, value)) = parse_xvg_metadata_line(&line_trimmed) {
                match key.as_str() {
                    "title" => metadata.title = Some(value),
                    "xaxis  label" => metadata.x_label = Some(value),
                    "yaxis  label" => metadata.y_label = Some(value),
                    _ => {}
                }
            }
        } else {
            let parts: Vec<&str> = line_trimmed.split_whitespace().collect();
            if parts.is_empty() {
                continue;
            }

            let parsed: Result<Vec<f64>, _> = parts.iter()
                .map(|s| s.parse::<f64>())
                .collect();

            if let Ok(values) = parsed {
                if data_columns.is_empty() {
                    for _ in 0..values.len() {
                        data_columns.push(Vec::new());
                    }
                }

                if values.len() != data_columns.len() {
                    eprintln!("Warning: inconsistent number of columns: line: {}", line_trimmed);
                    continue; 
                }
                for (i, val) in values.iter().enumerate() {
                    data_columns[i].push(*val);
                }
            }
        }
    }
    Ok(XvgData { metadata, data:data_columns })
}

fn main() {
    let argparses = App::new("Gromacs XVG to SVG Converter")
        .version("0.1")
        .author("gutao@westlake.edu.cn")
        .arg(
            Arg::new("input")
                .help("input xvg file path")
                .required(true)
                .index(1),
        )
        .arg(
            Arg::new("output")
                .help("output svg file path")
                .short('o')
                .long("output")
                .takes_value(true),
        )
        .arg(
            Arg::new("draw_columns")
                .help("Select draw columns")
                .short('n')
                .long("columns")
                .default_value("1")
                .takes_value(true),
        )
        .arg(
            Arg::new("legend")
                .help("legend name")
                .short('l')
                .long("legend")
                .default_value("RMSF")
                .takes_value(true),
        )
        .get_matches();
    
        
    let xvgpath=argparses.value_of("input").unwrap();
    let ncolumn = argparses.value_of("draw_columns").unwrap();
    let legend = argparses.value_of("legend").unwrap();
    use std::path::Path;
    let svgpath = argparses.value_of("output").map(String::from).unwrap_or_else(|| {
        let path = Path::new(xvgpath);
        let stem = match path.file_stem().and_then(|s| s.to_str()) {
            Some(stem_str) => stem_str,
            None => {
                panic!("Unable to extract stem from file path");
            }
        };
        format!("{}.svg", stem)
    });
    
    let xvg_data = read_xvg(xvgpath).unwrap();
    let x_label = xvg_data.metadata.x_label.unwrap_or("X label".to_string());
    let y_label = xvg_data.metadata.y_label.unwrap_or("Y label".to_string());
    let title = xvg_data.metadata.title.unwrap_or("Title".to_string());
    let data = xvg_data.data;
    let x_data=data[0].clone();
    let x_max= x_data.iter().cloned().max_by(|a, b| a.partial_cmp(b).unwrap());
    let x_min= x_data.iter().cloned().min_by(|a, b| a.partial_cmp(b).unwrap());
    let mut series = Vec::new();

    for i in 1..data.len() {
        let y_data=data[i].clone();
        let pairs: Vec<(f64, f64)> = x_data.iter()
            .zip(y_data.iter())
            .map(|(&x, &y)| (x, y))
            .collect();    
        series.push(
            Series::builder()
                .name(&legend)
                .color("Black")
                .data(pairs.clone())
                .line(Line::Solid)
                .build(),
        )
    };

    let plot = QuillPlot::builder()
        .dimensions((600, 400))
        .title(&title)
        .x_label(&x_label)
        .y_label(&y_label)
        .grid(Grid::None)
        .data(series)
        .build();
    
    plot.to_svg(&svgpath).unwrap();

    let idx: usize = ncolumn.parse().unwrap();
    let points: Vec<(f32, f32)> = x_data.iter()
        .zip(data[idx].iter())
        .map(|(&x, &y)| (x as f32, y as f32))   
        .collect(); 

    if let Some((Width(w), Height(h))) = terminal_size() {
        let width_val = 1.7f32 * w as f32;
        let height_val = 2.8f32 * h as f32;

        let width = (width_val.ceil() as u16).into();
        let height = (height_val.ceil() as u16).into();
        Chart::new(width,height,x_min.unwrap() as f32,x_max.unwrap() as f32)
            .lineplot(&Shape::Lines(&points))
            .display();
    } else {
        Chart::new(120,40,x_min.unwrap() as f32,x_max.unwrap() as f32)
            .lineplot(&Shape::Lines(&points))
            .display();
    };
}
