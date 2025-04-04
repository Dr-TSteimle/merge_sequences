// Copyright (C) 2025 [Thomas Steimlé]
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.

use std::ops::Range;

#[derive(Debug, Clone, PartialEq)]
pub enum Nucleotide {
    A,
    C,
    G,
    T,
}

#[derive(Debug, Clone, PartialEq)]
pub struct DNAString(pub Vec<Nucleotide>);

impl DNAString {
    pub fn new(input: Vec<u8>) -> DNAString {
        // replace by new from vec u8 or from string
        let mut seq: Vec<Nucleotide> = Vec::new();
        for nt in input {
            match nt {
                65u8 => seq.push(Nucleotide::A),
                67u8 => seq.push(Nucleotide::C),
                71u8 => seq.push(Nucleotide::G),
                84u8 => seq.push(Nucleotide::T),
                _ => panic!("Unknown nucleotid {}", String::from_utf8_lossy(&[nt])),
            }
        }
        DNAString(seq)
    }
    pub fn new_empty() -> DNAString {
        // replace by new
        DNAString(Vec::new())
    }
    pub fn as_string(&self) -> String {
        let mut str = String::new();

        for nt in self.0.iter() {
            match nt {
                Nucleotide::A => str.push('A'),
                Nucleotide::C => str.push('C'),
                Nucleotide::G => str.push('G'),
                Nucleotide::T => str.push('T'),
            }
        }
        str
    }
    pub fn as_string_2(&self) -> String {
        self.0
            .iter()
            .map(|nt| match nt {
                Nucleotide::A => "A",
                Nucleotide::C => "C",
                Nucleotide::G => "G",
                Nucleotide::T => "T",
            })
            .map(String::from)
            .collect::<Vec<String>>()
            .join("")
    }
    pub fn k_mers(&self, k: usize) -> Vec<DNAString> {
        let len = self.0.len();
        if k <= len {
            let mut res: Vec<DNAString> = Vec::new();
            let diff = (len - k) + 1;
            for i in 0..diff {
                res.push(DNAString(self.0[i..(i + k)].to_vec()));
            }
            res
        } else {
            panic!("{} read is too short {}", self.as_string(), len);
        }
    }
    pub fn compressed_k_mers(&self, k: usize) -> Vec<CompressedString> {
        let len = self.0.len();
        if k < len {
            let mut res: Vec<CompressedString> = Vec::new();
            let diff = (len - k) + 1;
            for i in 0..diff {
                res.push(CompressedString::new(DNAString(
                    self.0[i..(i + k)].to_vec(),
                )));
            }
            res
        } else {
            vec![CompressedString::new(self.clone())]
        }
    }

    pub fn merge(
        &mut self,
        input: &DNAString,
        kmer_len: usize,
        max_mismatches: i8,
        max_consecutive_mismatches: i8,
    ) -> Result<MergeResult, MergeErrors> {
        let input_kmers = &input.k_mers(kmer_len);
        // let max_consecutive_mismatches = 1i8;
        // let max_mismatches = 2i8;

        // let kmer_size = self.kmers[0].0.len();
        // let cluster_n_kmers = self.kmers.len();
        // let cluster_seq_len = self.kmers.len() + (kmer_size - 1);

        // let mut clust_ranges: Range<usize>;
        // let mut input_ranges: Range<usize>;

        // right overlap

        let mut kmers = self.k_mers(kmer_len);

        let res_opt = overlapping_ranges(
            &kmers,
            input_kmers,
            max_consecutive_mismatches,
            max_mismatches,
        );
        let res: Overlap;
        let mut input_on_sequence = true;

        // TODO : begin with smallest on biggest

        // Try input on sequence
        if res_opt.is_some() {
            res = res_opt.unwrap();
        } else {
            // try sequence on input
            let res_opt = overlapping_ranges(
                input_kmers,
                &kmers,
                max_consecutive_mismatches,
                max_mismatches,
            );
            if res_opt.is_some() {
                input_on_sequence = false;
                res = res_opt.unwrap(); // TODO : reverse the mismatches postions
                                        // take input as new sequence
                                        // self.0 = input.0.clone();
            } else {
                return Err(MergeErrors::NoMerge);
            }
        }

        if res.range_a.len() == input_kmers.len() {
            if input_on_sequence {
                // Input is inside sequence
                Ok(MergeResult::Inside(res))
            } else {
                // Input is outside sequence
                // take input as new sequence
                self.0 = input.0.clone();
                Ok(MergeResult::Outside(res))
            }
        } else {
            // Input is NOT inside sequence
            if res.range_b.start > 0 {
                // Add to the left
                if res.range_b.start > input_kmers.len() {
                    let seq = kmers_to_dnastring(input_kmers);
                    println!("{input_kmers:?}");
                    println!("{seq:?}");
                    println!("{res:?}");
                    return Err(MergeErrors::Unexpected);
                }
                let kmers_to_add_left = &input_kmers[0..res.range_b.start];

                kmers_to_add_left
                    .iter()
                    .rev()
                    .for_each(|it| kmers.insert(0, it.clone()));
                self.0 = kmers_to_dnastring(&kmers).0;
                return Ok(MergeResult::Left(res));
            }
            if res.range_b.end < input_kmers.len() {
                // println!("Add to the right");
                let kmers_to_add_right = &input_kmers[res.range_b.end..];

                kmers_to_add_right
                    .iter()
                    .for_each(|it| kmers.insert(kmers.len(), it.clone()));
                self.0 = kmers_to_dnastring(&kmers).0;
                return Ok(MergeResult::Right(res));
            }
            Err(MergeErrors::Unexpected)
        }
    }
}
#[derive(Debug, Clone, Copy)]
pub enum MergeErrors {
    NoMerge,
    Unexpected,
}
#[derive(Debug, Clone)]
pub enum MergeResult {
    Inside(Overlap),
    Outside(Overlap),
    Left(Overlap),
    Right(Overlap),
}
#[derive(PartialEq)]
pub struct CompressedString {
    seq: Vec<NucleotidFreq>,
}

impl CompressedString {
    pub fn new(seq: DNAString) -> CompressedString {
        let mut res: Vec<NucleotidFreq> = Vec::new();
        res.push(NucleotidFreq::new(seq.0[0].clone()));
        for i in 1..seq.0.len() {
            let nt = &seq.0[i];
            let mut last = res.pop().unwrap();
            // let last = &mut res[res.len()];
            if nt == &last.nucleotid {
                last.increment();
                res.push(last);
            } else {
                res.push(last);
                res.push(NucleotidFreq::new(nt.clone()));
            }
        }
        CompressedString { seq: res }
    }

    pub fn as_string(&self) -> String {
        self.seq
            .iter()
            .flat_map(|it| {
                let mut str = String::new();
                match it.nucleotid {
                    Nucleotide::A => str.push('A'),
                    Nucleotide::C => str.push('C'),
                    Nucleotide::G => str.push('G'),
                    Nucleotide::T => str.push('T'),
                };
                if it.freq > 1 {
                    vec![str, it.freq.to_string()]
                } else {
                    vec![str]
                }
            })
            .collect::<Vec<String>>()
            .join("")
    }
}

#[derive(PartialEq)]
struct NucleotidFreq {
    nucleotid: Nucleotide,
    freq: u8,
}
impl NucleotidFreq {
    fn new(nucleotid: Nucleotide) -> NucleotidFreq {
        NucleotidFreq { nucleotid, freq: 1 }
    }
    fn increment(&mut self) {
        self.freq += 1;
    }
}

fn kmers_to_dnastring(kmers: &[DNAString]) -> DNAString {
    let len = kmers.len();

    let mut dna_string = DNAString::new_empty();
    for (pos, kmer) in kmers.iter().enumerate() {
        if pos != len - 1 {
            dna_string.0.push(kmer.0[0].clone());
        } else {
            kmer.0.iter().for_each(|nt| dna_string.0.push(nt.clone()));
        }
    }
    dna_string
}
#[derive(Debug, Clone)]
pub struct Overlap {
    pub range_a: Range<usize>,
    pub range_b: Range<usize>,
    pub mismatches: Option<Vec<(usize, Nucleotide)>>,
}

pub fn overlapping_ranges(
    a_kmers: &[DNAString],
    b_kmers: &[DNAString],
    max_consecutive_mismatches: i8,
    max_mismatches: i8,
) -> Option<Overlap> {
    let kmer_size = b_kmers[0].0.len();

    // ensure both have same same kmers length
    assert_eq!(kmer_size, a_kmers[0].0.len());

    let a_full_len = a_kmers.len();
    let b_full_len = b_kmers.len();

    // merging b_kmers from the left of a_kmers
    // compare each non overlapping a_kmer to the first b_kmer

    let mut a_pos_match: Option<usize> = None;
    // let mut b_pos_match: Option<usize> = None;

    let mut a_range: Option<Range<usize>> = None;
    let mut b_range: Option<Range<usize>> = None;

    let mut mismatches_pos = None;
    'main: for (b_pos, b_kmer) in b_kmers.iter().enumerate() {
        if b_pos >= ((max_consecutive_mismatches as usize) + 1) * kmer_size {
            break;
        }

        // println!("Tryng to match b_kmer numero {}", b_pos);

        let b_seq = &b_kmer.0;

        'kmer_best_match: for (a_pos, a_kmer) in a_kmers.iter().enumerate() {
            let a_seq = &a_kmer.0;
            if a_seq == b_seq {
                // println!("Found a kmer perfect match at nt {}", a_pos);
                // println!("kmer : {}", a_kmer.as_string());
                a_pos_match = Some(a_pos);
                // b_pos_match = Some(b_pos);
                break 'kmer_best_match;
            }
        }

        // no perfect match (small reads)
        // if a_pos_match.is_none() {
        //     'kmer_match: for (a_pos, a_kmer) in a_kmers.iter().enumerate() {
        //         let a_seq = &a_kmer.0;
        //         let (test, _) = compare_with_tolerance(a_seq, b_seq, max_consecutive_mismatches, max_mismatches);

        //         if test {
        //             // println!("Found a kmer perfect match at nt {}", a_pos);
        //             // println!("kmer : {}", a_kmer.as_string());
        //             a_pos_match = Some(a_pos);
        //             // b_pos_match = Some(b_pos); // for clarity but to be removed idem b_pos
        //             break 'kmer_match ;
        //         }
        //     }
        // }

        // if found an a_kmer that match the b_kmer;
        if a_pos_match.is_some() {
            // should come back with tol
            // and all front should match with tolerance
            let a_pos_tmp = a_pos_match.unwrap();
            a_pos_match = None;

            // let a_dist_back = a_pos_tmp;
            // let b_dist_back = b_pos;

            let max_back_distance = if a_pos_tmp < b_pos { a_pos_tmp } else { b_pos };

            let iter_distance = max_back_distance;
            // println!("Possible backward steps {}", iter_distance);

            // loop {
            let a_iter_pos = a_pos_tmp - iter_distance;
            let b_iter_pos = b_pos - iter_distance;

            let a_dist_to_end = a_full_len - a_iter_pos;
            let b_dist_to_end = b_full_len - b_iter_pos;

            let overlap_len = if a_dist_to_end < b_dist_to_end {
                a_dist_to_end
            } else {
                b_dist_to_end
            };

            let a_iter_range = a_iter_pos..(a_iter_pos + overlap_len);
            let b_iter_range = b_iter_pos..(b_iter_pos + overlap_len);

            let a_seq_iter = kmers_to_dnastring(&a_kmers[a_iter_range.clone()]); // optimize by removing nt from previous DNAstr
            let b_seq_iter = kmers_to_dnastring(&b_kmers[b_iter_range.clone()]);

            let compare_res = compare_with_tolerance(
                &a_seq_iter.0,
                &b_seq_iter.0,
                max_consecutive_mismatches,
                max_mismatches,
            );
            // println!("compare_res {:?}", compare_res);

            if compare_res.0 {
                // println!("Match !");
                a_range = Some(a_iter_range);
                b_range = Some(b_iter_range);
                mismatches_pos = compare_res.1;
                break 'main;
            }
        }
    }

    a_range.zip(b_range).map(|(range_a, range_b)| Overlap {
        range_a,
        range_b,
        mismatches: mismatches_pos,
    })
}

// compare each nucleotids with sequence tolerance
// a is the reference
fn compare_with_tolerance(
    a_seq: &[Nucleotide],
    b_seq: &[Nucleotide],
    max_consecutive_mismatches: i8,
    max_mismatches: i8,
) -> (bool, Option<Vec<(usize, Nucleotide)>>) {
    assert!(
        max_consecutive_mismatches <= max_mismatches,
        "consecutive mismatches should be at max equal to maximum allowed mismatches"
    );
    let seq_len = a_seq.len();
    assert_eq!(seq_len, b_seq.len(), "sequence length should be equal");
    assert!(
        (max_mismatches as usize) <= seq_len,
        "sequence length should be superior to maximum allowed mismatches"
    );

    // init mismatches counters
    let mut n_mismatches = 0i8;
    let mut n_consecutive_mismatches = 0i8;
    let mut mismatches_pos: Vec<(usize, Nucleotide)> = Vec::new();

    let mut last_match = true;
    let mut test = true;
    for i in 0..seq_len {
        let curr_match = a_seq[i] == b_seq[i];
        if !curr_match {
            n_mismatches += 1;
            mismatches_pos.push((i, b_seq[i].clone()));
            if !last_match {
                n_consecutive_mismatches += 1;
            }
        }
        last_match = curr_match;

        test = n_mismatches <= max_mismatches
            && n_consecutive_mismatches <= max_consecutive_mismatches;
        if !test {
            break;
        }
    }
    let mismatches_pos_opt = if !mismatches_pos.is_empty() {
        Some(mismatches_pos)
    } else {
        None
    };
    (test, mismatches_pos_opt)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kmers() {
        let seq_a = "ATTTCACATCTGTAATACTCTGTCTCCTTGTTATAATTTCATTTACTAGTTT";
        let sequence_a = DNAString::new(seq_a.as_bytes().to_vec());
        let kmers = sequence_a.k_mers(50);

        assert_eq!(seq_a.len() - 50 + 1, kmers.len());
    }

    #[test]
    fn left() {
        let seq_a = "ATTTCACATCTGTAATACTCTGTCTCCTTGTTATAATTTCATTTACTAGTTATAATTTATAATGCAAACTGGATTGCAGCCCCAGTGCCAGGACTCAAATTATCCCAGAAATATAGGAAAAAGATCAACTCACGGGGCTCCACGAAGAGTT";
        let seq_b = "GGCCTATTTCACATCTGTAATACTCTGTCTCCTTGTTATAATTTCATTTACTAGTTATAATTTATAATGCAAACTGGATTGCAGCCCCAGTGCCAGGACTCAAATTATCCCAGAAATATAGGAAAAAGATCAACTCACGGGGCTCCACGAA";

        let mut sequence_a = DNAString::new(seq_a.as_bytes().to_vec());
        let sequence_b = DNAString::new(seq_b.as_bytes().to_vec());

        match sequence_a.merge(&sequence_b, 20, 2, 1) {
            Err(err) => panic!("{:?}", err),
            Ok(result) => match result {
                MergeResult::Left(_) => println!("{:?}", result),
                _ => panic!("Error {:?}", result),
            },
        }
        assert_eq!(
            DNAString::new("GGCCT".as_bytes().to_vec()).0,
            sequence_a
                .0
                .into_iter()
                .take(5)
                .collect::<Vec<Nucleotide>>()
        );
    }

    #[test]
    fn inside() {
        let kmer_len = 20;
        let max_mismatches = 2;
        let max_consecutive_mismatches = 1;

        let seq_a = "ACTCACGGGGCTCCACGAAGAGTTTGATGCCAGTTCCGAAAAACATCTGTCGGGTGTCCCGGCATCCCCCAGAGGAAAGGGATGATATACATCTTACATAAATTGCAGTCACAAGAAAGAAAGAAAAACAGAGATAGCCCAATAAATTAGGAAACCTTTACACTGTTTCGGCTTGATGAATTAGAAGACAGAATCCACTATAC";
        let mut sequence_a = DNAString::new(seq_a.as_bytes().to_vec());

        // no mismatch
        println!("Should be inside without mismatch -------------------------------------");
        let seq_c = "AAAACATCTGTCGGGTGTCCCGGCATCCCCCAGAGGAAAGGGATGATATACATCTTACATAAATTGCAGTCACAAGAAAGAAAGAAAAACAGAGATAGCCCAATAAATTAGGAAACCTTTACAC";
        let sequence_c = DNAString::new(seq_c.as_bytes().to_vec());
        let result = sequence_a.merge(
            &sequence_c,
            kmer_len,
            max_mismatches,
            max_consecutive_mismatches,
        );
        println!("{:?}", result);
        if let Ok(ref r) = result {
            match r {
                MergeResult::Inside(_) => (),
                _ => panic!("Error {:?}", result),
            }
        } else {
            panic!("Error {:?}", result)
        }

        // one mismatch on left
        println!("------------ Should be inside with one mismatch ------------------------");
        let seq_c = "TAAACATCTGTCGGGTGTCCCGGCATCCCCCAGAGGAAAGGGATGATATACATCTTACATAAATTGCAGTCACAAGAAAGAAAGAAAAACAGAGATAGCCCAATAAATTAGGAAACCTTTACAC";
        let sequence_c = DNAString::new(seq_c.as_bytes().to_vec());
        let result = sequence_a.merge(
            &sequence_c,
            kmer_len,
            max_mismatches,
            max_consecutive_mismatches,
        );
        println!("{:?}", result);
        if let Ok(ref r) = result {
            match r {
                MergeResult::Inside(ov) => {
                    if ov.mismatches.as_ref().unwrap().len() != 1 {
                        panic!("Error {:?}", result)
                    }
                }
                _ => panic!("Error {:?}", result),
            }
        } else {
            panic!("Error {:?}", result)
        }

        // // // 2 mismatches on left
        println!("------------ Should be inside with two mismatches ---------------------");
        let seq_c = "TTAACATCTGTCGGGTGTCCCGGCATCCCCCAGAGGAAAGGGATGATATACATCTTACATAAATTGCAGTCACAAGAAAGAAAGAAAAACAGAGATAGCCCAATAAATTAGGAAACCTTTACAC";
        let sequence_c = DNAString::new(seq_c.as_bytes().to_vec());
        let result = sequence_a.merge(
            &sequence_c,
            kmer_len,
            max_mismatches,
            max_consecutive_mismatches,
        );
        println!("{:?}", result);
        match result {
            Ok(ref r) => {
                if let MergeResult::Inside(ov) = r {
                    if ov.mismatches.as_ref().unwrap().len() != 2 {
                        panic!("Error {:?}", result)
                    }
                } else {
                    panic!("Error {:?}", result)
                }
            }
            _ => panic!("Error {:?}", result),
        }

        // // // > 2 mismatches on left
        println!("------------ > 2 mismatches -------------------------------------------");
        let seq_c = "TTTACATCTGTCGGGTGTCCCGGCATCCCCCAGAGGAAAGGGATGATATACATCTTACATAAATTGCAGTCACAAGAAAGAAAGAAAAACAGAGATAGCCCAATAAATTAGGAAACCTTTACAC";
        let sequence_c = DNAString::new(seq_c.as_bytes().to_vec());
        let result = sequence_a.merge(
            &sequence_c,
            kmer_len,
            max_mismatches,
            max_consecutive_mismatches,
        );
        println!("{:?}", result);
        if result.is_ok() {
            panic!("Error {:?}", result)
        }
    }
}
