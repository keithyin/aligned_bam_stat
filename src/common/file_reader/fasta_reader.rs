use std::{collections::HashMap, fs, io::{BufRead, BufReader}};

use crate::common::pb_tools::{get_spin_pb, DEFAULT_INTERVAL};


#[derive(Debug, PartialEq)]
enum LineTag {
    Header,
    Seq
}

pub struct FastaReaderIter<'a> {
    reader: &'a mut dyn BufRead,
    header: Option<String>
}

impl<'a> FastaReaderIter<'a>{
    pub fn new(reader: &'a mut dyn BufRead) -> Self {
        Self { reader, header: None}
    }

    fn read_one_line(&mut self) -> Option<(String, LineTag)> {

        let mut line = String::new();
        if let Ok(n) = self.reader.read_line(&mut line) {
            if n == 0 {
                return None;
            }
            line = line.trim().to_string();
            if  line.starts_with(">") {
                return Some((line, LineTag::Header));
            } else {
                return Some((line, LineTag::Seq));
            }
        } else {
            return None;
        }
    }

    pub fn get_qname2seq(self, show_pbar: bool) -> HashMap<String, String> {
        let mut result = HashMap::new();

        let mut pbar = None;
        if show_pbar {
            pbar = Some(get_spin_pb(format!(">> reading fasta"), DEFAULT_INTERVAL));
        }
        for item in self {
            let q_name = item.0.split_whitespace()
                .take(1)
                .map(|v| &v[1..])
                .collect::<Vec<_>>()[0]
                .to_string();
            result.insert(q_name, item.1);
            if let Some(ref pbar_) = pbar {
                pbar_.inc(1);
            }
        }

        if let Some(ref pbar_) = pbar {
            pbar_.finish();
        }

        result
    }
}

impl<'a> Iterator for FastaReaderIter<'a> {
    type Item = (String, String);
    fn next(&mut self) -> Option<Self::Item> {
        let mut seq_buffer = String::new();
        loop {

            if self.header.is_none() {
                if let Some((line, tag)) = self.read_one_line() {
                    if tag != LineTag::Header {
                        panic!("invalid");
                    }

                    self.header = Some(line);
                    continue;

                } else {
                    return None;
                }
            }

            if let Some((line, tag)) = self.read_one_line() {
                match tag {
                    LineTag::Seq => {seq_buffer.push_str(&line);},
                    LineTag::Header => {
                        let res = Some((self.header.take().unwrap(), seq_buffer));
                        self.header = Some(line);
                        return res;                        
                    }
                }
            } else {
                return Some((self.header.take().unwrap(), seq_buffer));
            }

        }
    }
}


pub struct FastaFile {
    ref_name2seq: HashMap<String, String>
}

impl FastaFile {
    
    pub fn new(filepath: &str) -> Self {
        let reader = fs::File::open(filepath).expect(&format!("filepath:'{}'", filepath));
        let mut buf_reader = BufReader::new(reader);
        let iter = FastaReaderIter::new(&mut buf_reader);
        Self { ref_name2seq: iter.get_qname2seq(false) }
    }

    pub fn get_ref_seq(&self, refname: &str) -> Option<&str> {
        self.ref_name2seq.get(refname).and_then(|v| Some(v.as_str()))
    }

    pub fn get_ref_name2seq(&self) -> &HashMap<String, String> {
        &self.ref_name2seq
    }

}


#[cfg(test)]
mod test {
    use std::{fs, io::BufReader};

    use super::FastaReaderIter;

    #[test]
    fn test_fasta_reader_iter() {
        let fasta_file = "/data/ccs_data/HG002/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta";
        let file = fs::File::open(fasta_file).unwrap();
        let mut reader = BufReader::new(file);

        let iter = FastaReaderIter::new(&mut reader);
        for (i, v) in iter.into_iter().enumerate() {
            eprintln!("{:?}", v.0);
            eprintln!("{:?}", &v.1.as_str()[..1000]);

        }

    }
}