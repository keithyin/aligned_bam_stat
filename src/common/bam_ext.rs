use rust_htslib::bam::{self, ext::BamRecordExtensions, record::Cigar, Record};
pub struct BamRecordExt<'a> {
    record: &'a Record
}

#[allow(unused)]
impl<'a> BamRecordExt<'a> {
    
    pub fn new(record: &'a Record) -> Self {
        Self { record }
    }

    pub fn get_qname(&self) -> String {
        String::from_utf8(self.record.qname().to_vec()).unwrap()
    }

    pub fn get_seq(&self) -> String {
        String::from_utf8(self.record.seq().as_bytes()).unwrap()
    }

    pub fn get_coverage(&self) -> Option<f32> {
        self.get_float(b"ec")
    }

    pub fn get_identity(&self) -> Option<f32> {
        self.get_float(b"iy")
    }

    pub fn get_rq(&self) -> Option<f32> {
        self.get_float(b"rq")

    }

    pub fn get_ch(&self) -> Option<usize> {
        self.get_int(b"ch")
    }


    pub fn get_np(&self) -> Option<usize> {
        self.get_int(b"np")
    }

    pub fn get_dw(&self) -> Option<Vec<usize>> {
        self.get_uint_list(b"dw")
    }

    pub fn get_cr(&self) -> Option<Vec<usize>> {
        self.get_uint_list(b"cr")
    }

    fn get_int(&self, tag: &[u8]) -> Option<usize> {
        self.record.aux(tag)
            .ok()
            .and_then(|aux| 
                match aux {
                    bam::record::Aux::I8(v) =>      Some(v as usize),
                    bam::record::Aux::U8(v) =>      Some(v as usize),
                    bam::record::Aux::I16(v) =>    Some(v as usize),
                    bam::record::Aux::U16(v) =>    Some(v as usize),
                    bam::record::Aux::I32(v) =>    Some(v as usize),
                    bam::record::Aux::U32(v) =>    Some(v as usize),
    
                    _ => None,
                }
            )
    }

    fn get_float(&self, tag: &[u8]) -> Option<f32> {
        self.record.aux(tag)
            .ok()
            .and_then(|aux| 
                match aux {
                    bam::record::Aux::Float(v) =>      Some(v as f32),
                    bam::record::Aux::Double(v) =>      Some(v as f32),
                    _ => None,
                }
            )
    }

    fn get_uint_list(&self, tag: &[u8]) -> Option<Vec<usize>> {
        self.record.aux(tag)
        .ok()
        .and_then(|aux| 
            match aux {
                bam::record::Aux::ArrayU8(v) => Some(v.iter().map(|v| v as usize).collect::<Vec<usize>>()),
                bam::record::Aux::ArrayU16(v) => Some(v.iter().map(|v| v as usize).collect::<Vec<usize>>()),
                bam::record::Aux::ArrayU32(v) => Some(v.iter().map(|v| v as usize).collect::<Vec<usize>>()),
                _ => None
            }
        )
    }


}


pub struct AlignedRecord<'a> {
    record: &'a Record
}

impl<'a> AlignedRecord<'a>{
    pub fn new(record: &'a Record) -> Self {
        Self { record }
    }

    /// matched / seq_len
    pub fn compute_effective_coverage(&self) -> f32 {
        let seq_len = self.record.seq_len_from_cigar(true);
        let matched = self.record.cigar().iter()
            .map(|cigar| {
                match cigar {
                    Cigar::Equal(n) | Cigar::Diff(n) | Cigar::Ins(n) => *n,
                    _ => 0     
                }
            })
            .sum::<u32>()
            ;

        if seq_len == 0 {
            0.0
        } else {
            matched as f32 / seq_len as f32
        }
    }

    pub fn compute_identity(&self) -> f32 {
        let tot_len = self.record.cigar().iter()
            .map(|cigar| {
                match cigar {
                    Cigar::Equal(n) | Cigar::Diff(n) | 
                    Cigar::Del(n) | Cigar::Ins(n) => *n,
                    
                    _ => 0
                }
            }).sum::<u32>();
        
        let matched = self.record.cigar().iter()
        .map(|cigar| {
            match cigar {
                Cigar::Equal(n) => *n,
                _ => 0
            }
        }).sum::<u32>();

        if tot_len == 0 {
            0.0
        } else {
            matched as f32 / tot_len as f32
        }

    }

}