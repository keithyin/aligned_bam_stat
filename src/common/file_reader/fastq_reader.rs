use std::io::BufRead;



// header, seq, qual
#[derive(Debug)]
pub struct FastqRecord(pub String, pub String, pub String);

impl FastqRecord {

    pub fn eqaul(&self, other: &FastqRecord) -> bool {
        self.1.eq(&other.1) && self.2.eq(&other.2)
    }

}

pub struct FastqReaderIter<'a> {
    reader: &'a mut dyn BufRead,
}

impl<'a> FastqReaderIter<'a>{
    pub fn new(reader: &'a mut dyn BufRead) -> Self {
        Self { reader}
    }

    fn read_one_line(&mut self) -> Option<String> {

        let mut line = String::new();
        if let Ok(n) = self.reader.read_line(&mut line) {
            if n == 0 {
                return None;
            }
            line = line.trim().to_string();
        } else {
            return None;
        }
        return Some(line);
    }

}

impl<'a> Iterator for FastqReaderIter<'a> {
    type Item = FastqRecord;
    fn next(&mut self) -> Option<Self::Item> {
        let header = self.read_one_line();
        if header.is_none() || header.as_ref().unwrap().trim().len() == 0{
            return None;
        }
        let mut header = header.unwrap();
        if !header.starts_with("@") {
            panic!("header:'{}' not a valid fastq header", header);
        }
        
        header = if let Some((v0, _)) = header.trim().split_once(" ") {
            v0[1..].to_string()
        } else {
            header[1..].to_string()
        };

        let seq = self.read_one_line().expect("not a valid FastqRecord").trim().to_string();
        let plus = self.read_one_line().expect("not a valid FastqRecord").trim().to_string();
        
        if !plus.starts_with("+") {
            panic!("plus:'{}' not a valid fastq plus", plus);
        }
        
        let qual = self.read_one_line().expect("not a valid FastqRecord").trim().to_string();

        Some(FastqRecord(header, seq, qual))

    }
}

#[cfg(test)]
mod test {
    use std::{fs, io::BufReader};

    use super::FastqReaderIter;

    #[test]
    fn test_fastq_reader_iter() {
        let fasta_file = "/data/ccs_data/ccs_eval_2024q1/icing_eval_20240506/20240329_Sync_Y0004_05_H01_Run0002_called_subreads.ccs_icing_result_staged/consensus.fastq";
        let file = fs::File::open(fasta_file).unwrap();
        let mut reader = BufReader::new(file);
        let iter = FastqReaderIter::new(&mut reader);
        for (i, v) in iter.into_iter().enumerate() {
            eprintln!("{:?}", i);
            eprintln!("{:?}", v);
            if i > 10 {
                break;
            }
        }

    }
}