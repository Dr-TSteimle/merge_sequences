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
//! # merge_sequences
//!
//! `merge_sequences` is a Rust library for approximate merging and alignment of DNA sequences,
//! using a fast k-mer based strategy with configurable mismatch tolerances.
//!
//! This library is ideal for applications in read clustering, contig extension, or
//! local alignment of noisy DNA sequences.
//!
//! ## 🧬 Key Types
//!
//! - [`DNAString`](struct.DNAString.html): Core DNA sequence container.
//! - [`Nucleotide`](enum.Nucleotide.html): Enum of valid DNA bases (A, C, G, T).
//! - [`MergeResult`](enum.MergeResult.html): Result of a merge (`Left`, `Right`, `Inside`, `Outside`).
//! - [`Overlap`](struct.Overlap.html): Describes the aligned regions and mismatches between sequences.
//!
//! ## ⚙️ Example
//!
//! ```rust
//! use merge_sequences::{DNAString, MergeResult};
//! use std::str::FromStr;
//!
//! let a = DNAString::from_str("ACGTACGT").unwrap();
//! let b = DNAString::from_str("GTACGAAA").unwrap();
//!
//! let mut merged = a.clone();
//! let result = merged.merge(&b, 4, 1, 1).unwrap();
//!
//! match result {
//!     MergeResult::Right(_) => {
//!         assert!(merged.as_string().starts_with("ACGT"));
//!         assert!(merged.as_string().ends_with("CGAAA"));
//!     },
//!     _ => panic!("Unexpected merge result: {:?}", result),
//! }
//! ```
//!
//! ## 📦 Installation
//!
//! Add this to your `Cargo.toml`:
//!
//! ```toml
//! [dependencies]
//! merge_sequences = "0.1"
//! ```
//!
//! ## © Thomas Steimlé
//!
//! Developed by [Thomas Steimlé](https://github.com/Dr-TSteimle)  
//! License: [GNU AGPL v3.0 or later](https://www.gnu.org/licenses/agpl-3.0.html)

use std::{cmp::min, ops::Range};

/// Represents a DNA nucleotide: A, C, G, or T.
///
/// This enum is used as the base unit in `DNAString` and provides conversion
/// from raw ASCII-encoded bytes (`u8`) and formatting support.
///
/// # Variants
///
/// - `A` - Adenine
/// - `C` - Cytosine
/// - `G` - Guanine
/// - `T` - Thymine
///
/// # Example
///
/// ```
/// use std::convert::TryFrom;
/// use merge_sequences::Nucleotide;
/// let nt = Nucleotide::try_from(b'G').unwrap();
/// assert_eq!(nt, Nucleotide::G);
///
/// let s = format!("{}", nt);
/// assert_eq!(s, "G");
/// ```
#[derive(Debug, Clone, PartialEq)]
pub enum Nucleotide {
    A,
    C,
    G,
    T,
}

impl TryFrom<u8> for Nucleotide {
    type Error = String;

    /// Converts a `u8` (ASCII-encoded) into a `Nucleotide`.
    ///
    /// # Errors
    ///
    /// Returns an error if the byte does not represent one of `A`, `C`, `G`, or `T`.
    fn try_from(nt: u8) -> Result<Self, Self::Error> {
        match nt {
            b'A' => Ok(Nucleotide::A),
            b'C' => Ok(Nucleotide::C),
            b'G' => Ok(Nucleotide::G),
            b'T' => Ok(Nucleotide::T),
            _ => Err(format!("Unknown nucleotide: {}", nt as char)),
        }
    }
}

impl std::fmt::Display for Nucleotide {
    /// Formats the `Nucleotide` as a single character (`A`, `C`, `G`, `T`).
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Nucleotide::A => 'A',
                Nucleotide::C => 'C',
                Nucleotide::G => 'G',
                Nucleotide::T => 'T',
            }
        )
    }
}

/// A wrapper around a vector of `Nucleotide`, representing a DNA sequence.
///
/// `DNAString` is the fundamental sequence container for DNA-based data structures.
/// It supports construction from a string, from a vector of ASCII bytes, and from iterators of `Nucleotide`.
///
/// # Examples
///
/// ## From string:
/// ```
/// use std::str::FromStr;
/// use merge_sequences::{DNAString, Nucleotide};
///
/// let dna = DNAString::from_str("ACGT").unwrap();
/// assert_eq!(dna.0.len(), 4);
/// ```
///
/// ## From byte vector:
/// ```
/// use std::convert::TryFrom;
/// use merge_sequences::{DNAString, Nucleotide};
///
/// let dna = DNAString::try_from(vec![b'A', b'C', b'G', b'T']).unwrap();
/// assert_eq!(dna.0[2], Nucleotide::G);
/// ```
///
/// ## From iterator:
/// ```
/// use merge_sequences::{DNAString, Nucleotide};
///
/// let nts = vec![Nucleotide::A, Nucleotide::C, Nucleotide::T];
/// let dna: DNAString = nts.into_iter().collect();
/// assert_eq!(dna.0.len(), 3);
/// ```
#[derive(Debug, Clone, PartialEq, Default)]
pub struct DNAString(pub Vec<Nucleotide>);

impl FromIterator<Nucleotide> for DNAString {
    fn from_iter<I: IntoIterator<Item = Nucleotide>>(iter: I) -> Self {
        DNAString(iter.into_iter().collect())
    }
}

impl TryFrom<Vec<u8>> for DNAString {
    type Error = String;

    fn try_from(input: Vec<u8>) -> Result<Self, Self::Error> {
        input.into_iter().map(Nucleotide::try_from).collect()
    }
}

impl std::str::FromStr for DNAString {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        s.bytes().map(Nucleotide::try_from).collect()
    }
}

impl DNAString {
    /// Constructs a new `DNAString` from a vector of ASCII-encoded nucleotides (`A`, `C`, `G`, `T`).
    ///
    /// # Panics
    ///
    /// Panics if the input contains invalid nucleotide codes.
    ///
    /// # Examples
    ///
    /// ```
    /// use merge_sequences::{ Nucleotide, DNAString };
    ///
    /// let dna = DNAString::new(b"ACGT".to_vec());
    /// ```
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

    /// Converts the `DNAString` back into a `String` of nucleotide characters.
    ///
    /// This is the inverse of constructing a `DNAString` from a `&str` or `Vec<u8>`.
    ///
    /// # Example
    /// ```
    /// use merge_sequences::DNAString;
    /// use std::str::FromStr;
    ///
    /// let dna = DNAString::from_str("ACGT").unwrap();
    /// assert_eq!(dna.as_string(), "ACGT");
    /// ```
    pub fn as_string(&self) -> String {
        self.0.iter().map(ToString::to_string).collect()
    }

    /// Splits the sequence into all overlapping k-mers of length `k`.
    ///
    /// The resulting `Vec<DNAString>` contains all substrings of length `k`
    /// that span the sequence, sliding one base at a time.
    ///
    /// # Panics
    /// Panics if `k > self.len()`
    ///
    /// # Example
    /// ```
    /// use merge_sequences::DNAString;
    /// use std::str::FromStr;
    ///
    /// let dna = DNAString::from_str("ACGT").unwrap();
    /// let kmers = dna.k_mers(2);
    /// let kmers_str: Vec<String> = kmers.iter().map(|k| k.as_string()).collect();
    /// assert_eq!(kmers_str, vec!["AC", "CG", "GT"]);
    /// ```
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

    /// Returns an iterator over all k-mers of length `k` as `DNAString`s.
    ///
    /// This is more efficient than `k_mers()` when the result does not need to be stored in memory.
    ///
    /// # Example
    /// ```
    /// use merge_sequences::DNAString;
    /// use std::str::FromStr;
    ///
    /// let dna = DNAString::from_str("ACGT").unwrap();
    /// let mut iter = dna.k_mer_iter(2);
    /// assert_eq!(iter.next().unwrap().as_string(), "AC");
    /// assert_eq!(iter.next().unwrap().as_string(), "CG");
    /// assert_eq!(iter.next().unwrap().as_string(), "GT");
    /// assert!(iter.next().is_none());
    /// ```
    pub fn k_mer_iter(&self, k: usize) -> impl Iterator<Item = DNAString> + '_ {
        self.0.windows(k).map(|win| DNAString(win.to_vec()))
    }

    pub fn k_mers_vec(&self, k: usize) -> Vec<DNAString> {
        self.k_mer_iter(k).collect()
    }

    /// Attempts to merge the current `DNAString` with another (`input`) based on approximate k-mer overlaps.
    ///
    /// The method identifies the best overlapping region between the two sequences using shared k-mers,
    /// allowing for a configurable number of mismatches and consecutive mismatches. Depending on the type of
    /// overlap found, it either incorporates the input into the current sequence or replaces it entirely.
    ///
    /// # Arguments
    ///
    /// * `input` - The `DNAString` to merge into `self`.
    /// * `kmer_len` - Length of k-mers used for overlap detection.
    /// * `max_mismatches` - Maximum total mismatches allowed in the overlapping region.
    /// * `max_consecutive_mismatches` - Maximum allowed consecutive mismatches.
    ///
    /// # Returns
    ///
    /// A `Result` containing the merge outcome:
    /// - `Ok(MergeResult)` on success
    /// - `Err(MergeErrors)` if no acceptable overlap is found or an unexpected case is encountered.
    ///
    /// # Merge outcomes
    ///
    /// - `MergeResult::Inside`: The input is fully contained within `self`, so no extension is needed.
    /// - `MergeResult::Outside`: The input contains `self`, so `self` is replaced with `input`.
    /// - `MergeResult::Left`: The input overlaps on the left of `self` and is prepended.
    /// - `MergeResult::Right`: The input overlaps on the right of `self` and is appended.
    ///
    /// # Errors
    ///
    /// - `MergeErrors::NoMerge`: No valid overlapping region was found.
    /// - `MergeErrors::Unexpected`: Unexpected alignment geometry encountered (e.g., invalid overlap range).
    ///
    /// # Panics
    ///
    /// This method does not panic.
    ///
    /// # Example
    ///
    /// ```
    /// use std::str::FromStr;
    /// use merge_sequences::DNAString;
    ///
    /// let mut seq1 = DNAString::from_str("ACGTAC").unwrap();
    /// let seq2 = DNAString::from_str("TACGGA").unwrap();
    /// let result = seq1.merge(&seq2, 3, 1, 1);
    /// ```
    pub fn merge(
        &mut self,
        input: &DNAString,
        kmer_len: usize,
        max_mismatches: i8,
        max_consecutive_mismatches: i8,
    ) -> Result<MergeResult, MergeErrors> {
        let self_kmers = self.k_mers_vec(kmer_len);
        let input_kmers = input.k_mers_vec(kmer_len);

        let forward = overlapping_ranges(
            &self_kmers,
            &input_kmers,
            max_consecutive_mismatches,
            max_mismatches,
        );
        let reverse = overlapping_ranges(
            &input_kmers,
            &self_kmers,
            max_consecutive_mismatches,
            max_mismatches,
        );

        let (overlap, input_on_sequence) = match (forward, reverse) {
            (Some(r1), Some(r2)) => {
                if r1.range_a.len() >= r2.range_b.len() {
                    (r1, true)
                } else {
                    (r2, false)
                }
            }
            (Some(r1), None) => (r1, true),
            (None, Some(r2)) => (r2, false),
            (None, None) => return Err(MergeErrors::NoMerge),
        };

        let input_len = input_kmers.len();
        let self_len = self_kmers.len();

        // ✅ Force containment replacement even if overlap was detected in forward direction
        if overlap.range_b.len() == self_len && input.0.len() > self.0.len() {
            self.0 = input.0.clone();
            return Ok(MergeResult::Outside(overlap));
        }

        if overlap.range_a.len() == input_len && input_on_sequence {
            return Ok(MergeResult::Inside(overlap));
        }

        let mut merged_kmers = self_kmers;

        // Extend to the left
        if overlap.range_b.start > 0 {
            let extra_left = &input_kmers[0..overlap.range_b.start];
            for kmer in extra_left.iter().rev() {
                merged_kmers.insert(0, kmer.clone());
            }

            self.0 = kmers_to_dnastring(&merged_kmers).0;
            return Ok(MergeResult::Left(overlap));
        }

        // Extend to the right
        if overlap.range_b.end < input_len {
            let extra_right = &input_kmers[overlap.range_b.end..];
            merged_kmers.extend(extra_right.iter().cloned());

            self.0 = kmers_to_dnastring(&merged_kmers).0;
            return Ok(MergeResult::Right(overlap));
        }

        Err(MergeErrors::Unexpected)
    }

    pub fn merge_old(
        &mut self,
        input: &DNAString,
        kmer_len: usize,
        max_mismatches: i8,
        max_consecutive_mismatches: i8,
    ) -> Result<MergeResult, MergeErrors> {
        let input_kmers = input.k_mers(kmer_len);
        let self_kmers = self.k_mers(kmer_len);

        // Try self vs input
        let mut input_on_sequence = true;
        let mut overlap = overlapping_ranges(
            &self_kmers,
            &input_kmers,
            max_consecutive_mismatches,
            max_mismatches,
        );

        // Try reversed direction if no overlap
        if overlap.is_none() {
            overlap = overlapping_ranges(
                &input_kmers,
                &self_kmers,
                max_consecutive_mismatches,
                max_mismatches,
            );
            if overlap.is_some() {
                input_on_sequence = false;
            } else {
                return Err(MergeErrors::NoMerge);
            }
        }

        let overlap = overlap.unwrap();

        if overlap.range_a.len() == input_kmers.len() {
            if input_on_sequence {
                Ok(MergeResult::Inside(overlap))
            } else {
                self.0 = input.0.clone();
                Ok(MergeResult::Outside(overlap))
            }
        } else {
            let mut kmers = self_kmers;

            // Extend left
            if overlap.range_b.start > 0 {
                if overlap.range_b.start > input_kmers.len() {
                    return Err(MergeErrors::Unexpected);
                }

                let left_extension = &input_kmers[..overlap.range_b.start];
                for kmer in left_extension.iter().rev() {
                    kmers.insert(0, kmer.clone());
                }

                self.0 = kmers_to_dnastring(&kmers).0;
                return Ok(MergeResult::Left(overlap));
            }

            // Extend right
            if overlap.range_b.end < input_kmers.len() {
                let right_extension = &input_kmers[overlap.range_b.end..];
                kmers.extend(right_extension.iter().cloned());

                self.0 = kmers_to_dnastring(&kmers).0;
                return Ok(MergeResult::Right(overlap));
            }

            Err(MergeErrors::Unexpected)
        }
    }
}

/// Represents possible errors that can occur during a DNA sequence merge operation.
///
/// These errors are returned by the [`DNAString::merge`] method when merging fails
/// due to non-overlapping sequences or unexpected edge cases.
///
/// # Variants
///
/// - `NoMerge`: No valid overlap was found between the sequences.
/// - `Unexpected`: An unexpected or inconsistent overlap configuration was detected
///   (e.g., mismatch between computed ranges and k-mer assumptions).
///
/// # Example
/// ```
/// use merge_sequences::{DNAString, MergeErrors};
/// use std::str::FromStr;
///
/// let mut a = DNAString::from_str("AAAAAA").unwrap();
/// let b = DNAString::from_str("TTTTTT").unwrap();
///
/// let result = a.merge(&b, 3, 1, 1);
/// assert!(matches!(result, Err(MergeErrors::NoMerge)));
/// ```
#[derive(Debug, Clone, Copy)]
pub enum MergeErrors {
    /// No valid overlap could be found.
    NoMerge,
    /// An internal inconsistency occurred during merging.
    Unexpected,
}

/// Represents the result of a successful merge between two `DNAString`s.
///
/// Returned by the [`DNAString::merge`] method, this enum indicates how
/// the merge was performed:
///
/// - `Inside`: The input was already fully contained within the existing sequence.
/// - `Outside`: The existing sequence was fully contained in the input, and replaced by it.
/// - `Left`: The input extended the sequence by overlapping on the left.
/// - `Right`: The input extended the sequence by overlapping on the right.
///
/// Each variant carries an [`Overlap`] struct with the matching ranges and optional mismatch data.
///
/// # Example
/// ```
/// use merge_sequences::{DNAString, MergeResult};
/// use std::str::FromStr;
///
/// let a = DNAString::from_str("ACGTACGT").unwrap();
/// let b = DNAString::from_str("GTACGAAA").unwrap();
///
/// let mut merged = a.clone();
/// let result = merged.merge(&b, 4, 1, 1).unwrap();
///
/// match result {
///     MergeResult::Right(overlap) => {
///         assert!(merged.as_string().contains("ACGTACGAAA"));
///     }
///     _ => panic!("Unexpected result: {:?}", result),
/// }
/// ```
#[derive(Debug, Clone)]
pub enum MergeResult {
    /// The input was fully contained within the existing sequence.
    Inside(Overlap),
    /// The existing sequence was fully contained in the input; it was replaced.
    Outside(Overlap),
    /// The input extended the sequence on the left side.
    Left(Overlap),
    /// The input extended the sequence on the right side.
    Right(Overlap),
}

// #[derive(PartialEq)]
// pub struct CompressedString {
//     seq: Vec<NucleotidFreq>,
// }
//
// impl CompressedString {
//     pub fn new(seq: DNAString) -> CompressedString {
//         let mut res: Vec<NucleotidFreq> = Vec::new();
//         res.push(NucleotidFreq::new(seq.0[0].clone()));
//         for i in 1..seq.0.len() {
//             let nt = &seq.0[i];
//             let mut last = res.pop().unwrap();
//             // let last = &mut res[res.len()];
//             if nt == &last.nucleotid {
//                 last.increment();
//                 res.push(last);
//             } else {
//                 res.push(last);
//                 res.push(NucleotidFreq::new(nt.clone()));
//             }
//         }
//         CompressedString { seq: res }
//     }
//
//     pub fn as_string(&self) -> String {
//         self.seq
//             .iter()
//             .flat_map(|it| {
//                 let mut str = String::new();
//                 match it.nucleotid {
//                     Nucleotide::A => str.push('A'),
//                     Nucleotide::C => str.push('C'),
//                     Nucleotide::G => str.push('G'),
//                     Nucleotide::T => str.push('T'),
//                 };
//                 if it.freq > 1 {
//                     vec![str, it.freq.to_string()]
//                 } else {
//                     vec![str]
//                 }
//             })
//             .collect::<Vec<String>>()
//             .join("")
//     }
// }

/// Represents the frequency count of a nucleotide in a collection.
///
/// This structure is useful when computing consensus sequences or tracking base support
/// across aligned reads or overlapping k-mers.
#[derive(PartialEq)]
pub struct NucleotidFreq {
    /// The nucleotide being counted.
    pub nucleotid: Nucleotide,
    /// The number of occurrences observed.
    pub freq: u8,
}

impl NucleotidFreq {
    /// Creates a new `NucleotidFreq` with frequency initialized to 1.
    ///
    /// # Arguments
    ///
    /// * `nucleotid` - The nucleotide to be counted.
    ///
    /// # Example
    ///
    /// ```
    /// use merge_sequences::{NucleotidFreq, Nucleotide};
    ///
    /// let freq = NucleotidFreq::new(Nucleotide::A);
    /// assert_eq!(freq.freq, 1);
    /// ```
    pub fn new(nucleotid: Nucleotide) -> NucleotidFreq {
        NucleotidFreq { nucleotid, freq: 1 }
    }

    /// Increments the frequency count by 1.
    ///
    /// # Example
    ///
    /// ```
    /// use merge_sequences::{NucleotidFreq, Nucleotide};
    ///
    /// let mut freq = NucleotidFreq::new(Nucleotide::C);
    /// freq.increment();
    /// assert_eq!(freq.freq, 2);
    /// ```
    pub fn increment(&mut self) {
        self.freq += 1;
    }
}

/// Reconstructs a full `DNAString` sequence from a list of overlapping k-mers.
///
/// This function assumes that the k-mers are consecutive and overlapping by `k-1` bases,
/// as produced by the `DNAString::k_mers()` method. The result is a single contiguous sequence.
///
/// # Arguments
///
/// * `kmers` - A slice of overlapping `DNAString` k-mers (e.g., of length `k`).
///
/// # Returns
///
/// A new `DNAString` representing the reconstructed sequence.
///
/// # Example
///
/// ```
/// use std::str::FromStr;
/// use merge_sequences::{DNAString, kmers_to_dnastring};
///
/// let kmers = vec![
///     DNAString::from_str("ACG").unwrap(),
///     DNAString::from_str("CGT").unwrap(),
///     DNAString::from_str("GTT").unwrap(),
/// ];
/// let sequence = kmers_to_dnastring(&kmers);
/// assert_eq!(sequence.as_string(), "ACGTT");
/// ```
///
/// # Panics
///
/// Panics if the input slice is empty or contains non-overlapping k-mers.
pub fn kmers_to_dnastring(kmers: &[DNAString]) -> DNAString {
    let len = kmers.len();

    let mut dna_string = DNAString::default();
    for (pos, kmer) in kmers.iter().enumerate() {
        if pos != len - 1 {
            dna_string.0.push(kmer.0[0].clone());
        } else {
            kmer.0.iter().for_each(|nt| dna_string.0.push(nt.clone()));
        }
    }
    dna_string
}

/// Represents an overlapping region between two sequences identified by shared or similar k-mers.
///
/// This structure is typically returned by functions like `overlapping_ranges()` and
/// describes the aligned regions (in k-mer space) between two sequences `A` and `B`.
///
/// It also optionally records mismatch positions if approximate matching was performed.
///
/// # Fields
///
/// - `range_a`: The inclusive range of overlapping k-mers in sequence A.
/// - `range_b`: The inclusive range of overlapping k-mers in sequence B.
/// - `mismatches`: Optional vector of mismatch positions and expected nucleotides
///   relative to the overlap alignment.
///
/// # Example
///
/// ```
/// use merge_sequences::{Overlap, Nucleotide};
///
/// let overlap = Overlap {
///     range_a: 2..5,
///     range_b: 0..3,
///     mismatches: Some(vec![(4, Nucleotide::T)]),
/// };
/// assert_eq!(overlap.range_a.len(), overlap.range_b.len());
/// ```
#[derive(Debug, Clone)]
pub struct Overlap {
    /// Range of overlapping k-mers in the first sequence (`a_kmers`)
    pub range_a: Range<usize>,
    /// Range of overlapping k-mers in the second sequence (`b_kmers`)
    pub range_b: Range<usize>,
    /// Optional list of mismatch positions and expected nucleotides
    pub mismatches: Option<Vec<(usize, Nucleotide)>>,
}

/// Searches for an approximate overlapping region between two DNA sequences,
/// given as slices of k-mers, using a seed-and-extend strategy.
///
/// This function identifies a seed k-mer in `b_kmers` that exactly matches one in `a_kmers`,
/// then attempts to extend that region in both directions (backwards and forwards)
/// to validate a full overlap under mismatch constraints.
///
/// The comparison tolerates a configurable number of total and consecutive mismatches,
/// as determined by `compare_with_tolerance()`.
///
/// # Arguments
///
/// * `a_kmers` - The k-merized form of the first sequence (usually the existing reference).
/// * `b_kmers` - The k-merized form of the second sequence (usually the new candidate to merge).
/// * `max_consecutive_mismatches` - Maximum number of allowed consecutive mismatches in the overlap.
/// * `max_mismatches` - Maximum total number of mismatches allowed in the overlap.
///
/// # Returns
///
/// An `Option<Overlap>` indicating the matched regions if a valid overlap is found.
/// The returned `Overlap` contains:
/// - `range_a`: The overlapping range in `a_kmers`.
/// - `range_b`: The overlapping range in `b_kmers`.
/// - `mismatches`: Optional positions of mismatches (if applicable).
///
/// Returns `None` if no valid overlap is found.
///
/// # Panics
///
/// Panics if the first k-mer in `a_kmers` or `b_kmers` has a different nucleotide length.
///
/// # Assumptions
///
/// * Both inputs are non-empty and contain k-mers of the same length.
/// * `kmers_to_dnastring()` and `compare_with_tolerance()` are consistent with
///   how the original k-mers were created.
///
/// # Example
///
/// ```
/// use std::str::FromStr;
/// use merge_sequences::{DNAString, overlapping_ranges};
///
/// // Helper to generate k-mers
/// fn kmerize(seq: &str, k: usize) -> Vec<DNAString> {
///     DNAString::from_str(seq).unwrap().k_mers(k)
/// }
///
/// let a = kmerize("AACGT", 3);  // ["AAC", "ACG", "CGT"]
/// let b = kmerize("CGTTT", 3);  // ["CGT", "GTT", "TTT"]
///
/// let overlap = overlapping_ranges(&a, &b, 1, 2);
///
/// assert!(overlap.is_some());
/// let overlap = overlap.unwrap();
///
/// assert_eq!(overlap.range_a, 2..3);  // "CGT"
/// assert_eq!(overlap.range_b, 0..1);  // "CGT"
/// ```
pub fn overlapping_ranges(
    a_kmers: &[DNAString],
    b_kmers: &[DNAString],
    max_consecutive_mismatches: i8,
    max_mismatches: i8,
) -> Option<Overlap> {
    let kmer_size = b_kmers[0].0.len();
    assert_eq!(kmer_size, a_kmers[0].0.len());

    let a_len = a_kmers.len();
    let b_len = b_kmers.len();

    let mut a_range: Option<Range<usize>> = None;
    let mut b_range: Option<Range<usize>> = None;
    let mut mismatch_positions = None;

    // For each b_kmer (up to a distance threshold), try to find a perfect match in a_kmers
    for (b_pos, b_kmer) in b_kmers.iter().enumerate() {
        if b_pos >= ((max_consecutive_mismatches as usize) + 1) * kmer_size {
            break;
        }

        if let Some(a_pos) = a_kmers.iter().position(|a_kmer| a_kmer.0 == b_kmer.0) {
            // Found a seed match — align backward to find maximum overlap window
            let max_back = min(a_pos, b_pos);
            let a_start = a_pos - max_back;
            let b_start = b_pos - max_back;

            let a_remaining = a_len - a_start;
            let b_remaining = b_len - b_start;
            let overlap_len = min(a_remaining, b_remaining);

            let a_range_tmp = a_start..(a_start + overlap_len);
            let b_range_tmp = b_start..(b_start + overlap_len);

            let a_seq = kmers_to_dnastring(&a_kmers[a_range_tmp.clone()]);
            let b_seq = kmers_to_dnastring(&b_kmers[b_range_tmp.clone()]);

            let (ok, mismatches) = compare_with_tolerance(
                &a_seq.0,
                &b_seq.0,
                max_consecutive_mismatches,
                max_mismatches,
            );

            if ok {
                a_range = Some(a_range_tmp);
                b_range = Some(b_range_tmp);
                mismatch_positions = mismatches;
                break;
            }
        }
    }

    a_range.zip(b_range).map(|(a, b)| Overlap {
        range_a: a,
        range_b: b,
        mismatches: mismatch_positions,
    })
}

/// Compares two nucleotide sequences with mismatch tolerance, enforcing both total and consecutive mismatch limits.
///
/// This function performs a base-by-base comparison between two equal-length slices of `Nucleotide`.
/// It returns whether the sequences are similar enough according to:
///
/// - `max_mismatches`: the maximum number of total mismatches allowed.
/// - `max_consecutive_mismatches`: the number of *consecutive* mismatches allowed,
///    but note: this counter is increased only when **two or more mismatches occur in a row**.
///
/// ### Key behavior:
///
/// - The **first mismatch in a run** does **not** increase the `n_consecutive_mismatches` counter.
/// - The counter is incremented only when **a mismatch follows another mismatch**.
/// - Example:
///     - mismatch → `n_consec = 0`
///     - mismatch, mismatch → `n_consec = 1`
///     - mismatch, mismatch, mismatch → `n_consec = 2`
///
/// This logic allows isolated mismatches without penalty, and starts counting consecutive mismatches
/// only from the second in a run. This is a crucial distinction from naive run-length counting.
///
/// # Arguments
///
/// * `a_seq` – Reference sequence, as a slice of `Nucleotide`.
/// * `b_seq` – Query sequence, same length as `a_seq`.
/// * `max_consecutive_mismatches` – Maximum number of *consecutive* mismatches (starting from the second in a run).
/// * `max_mismatches` – Maximum number of total mismatches allowed.
///
/// # Returns
///
/// Returns a tuple:
///
/// - `bool`: `true` if the sequences match within tolerance; `false` otherwise.
/// - `Option<Vec<(usize, Nucleotide)>>`: If mismatches were found (but tolerated), this contains the
///   0-based positions and expected `Nucleotide` values from `b_seq`. Returns `None` if sequences match perfectly.
///
/// # Panics
///
/// This function will panic if:
/// - `a_seq` and `b_seq` have different lengths.
/// - `max_consecutive_mismatches > max_mismatches`.
/// - `max_mismatches > a_seq.len()`.
///
/// # Example
///
/// ```
/// use merge_sequences::{Nucleotide, compare_with_tolerance};
///
/// let a = vec![Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T];
/// let b = vec![Nucleotide::A, Nucleotide::G, Nucleotide::G, Nucleotide::T]; // mismatch at position 1
///
/// let (ok, mismatches) = compare_with_tolerance(&a, &b, 1, 2);
/// assert!(ok);
/// assert_eq!(mismatches, Some(vec![(1, Nucleotide::G)]));
///
/// let (ok2, _) = compare_with_tolerance(&a, &b, 0, 0);
/// assert!(!ok2);
/// ```
pub fn compare_with_tolerance(
    a_seq: &[Nucleotide],
    b_seq: &[Nucleotide],
    max_consecutive_mismatches: i8,
    max_mismatches: i8,
) -> (bool, Option<Vec<(usize, Nucleotide)>>) {
    assert!(
        max_consecutive_mismatches <= max_mismatches,
        "consecutive mismatches should be at most equal to the maximum allowed mismatches"
    );
    assert_eq!(a_seq.len(), b_seq.len(), "sequence lengths must be equal");
    assert!(
        (max_mismatches as usize) <= a_seq.len(),
        "sequence length must be at least the number of allowed mismatches"
    );

    let seq_len = a_seq.len();
    let mut n_mismatches = 0i8;
    let mut n_consecutive_mismatches = 0i8;
    let mut mismatches_pos: Vec<(usize, Nucleotide)> = Vec::new();

    let mut last_match = true;

    for i in 0..seq_len {
        let curr_match = a_seq[i] == b_seq[i];

        if !curr_match {
            n_mismatches += 1;
            mismatches_pos.push((i, b_seq[i].clone()));

            if !last_match {
                n_consecutive_mismatches += 1;
            } else {
                n_consecutive_mismatches = 0; // ← important: only count consecutive mismatches beyond the first
            }
        }

        last_match = curr_match;

        if n_mismatches > max_mismatches || n_consecutive_mismatches > max_consecutive_mismatches {
            return (false, None);
        }
    }

    let mismatches_pos_opt = if mismatches_pos.is_empty() {
        None
    } else {
        Some(mismatches_pos)
    };

    (true, mismatches_pos_opt)
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
    fn perfect_match() {
        let a = vec![Nucleotide::A, Nucleotide::C, Nucleotide::G];
        let b = vec![Nucleotide::A, Nucleotide::C, Nucleotide::G];

        let (ok, mismatches) = compare_with_tolerance(&a, &b, 1, 2);
        assert!(ok);
        assert!(mismatches.is_none());
    }

    #[test]
    fn single_mismatch_allowed() {
        let a = vec![Nucleotide::A, Nucleotide::C, Nucleotide::G];
        let b = vec![Nucleotide::A, Nucleotide::T, Nucleotide::G];

        let (ok, mismatches) = compare_with_tolerance(&a, &b, 1, 1);
        assert!(ok);
        assert_eq!(mismatches, Some(vec![(1, Nucleotide::T)]));
    }

    #[test]
    fn two_consecutive_mismatches_blocked() {
        let a = vec![Nucleotide::A, Nucleotide::C, Nucleotide::G];
        let b = vec![Nucleotide::A, Nucleotide::T, Nucleotide::T];

        let (ok, _) = compare_with_tolerance(&a, &b, 0, 2);
        assert!(!ok); // two mismatches in a row, but only 0 consecutive allowed
    }

    #[test]
    fn total_mismatch_exceeds_limit() {
        let a = vec![Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T];
        let b = vec![Nucleotide::T, Nucleotide::T, Nucleotide::T, Nucleotide::T];

        let (ok, _) = compare_with_tolerance(&a, &b, 2, 2);
        assert!(!ok); // 3 mismatches > max allowed
    }

    #[test]
    fn consecutive_logic_exact() {
        let a = vec![Nucleotide::A, Nucleotide::T, Nucleotide::T, Nucleotide::G];
        let b = vec![Nucleotide::C, Nucleotide::C, Nucleotide::T, Nucleotide::G];

        // mismatch at 0 (isolated), mismatch at 1 follows a mismatch → n_consec = 1
        let (ok, _) = compare_with_tolerance(&a, &b, 1, 3);
        assert!(ok);

        let (fail, _) = compare_with_tolerance(&a, &b, 0, 3);
        assert!(!fail);
    }

    #[test]
    fn two_mismatches_separated_allowed() {
        let a = vec![Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T];
        let b = vec![Nucleotide::T, Nucleotide::C, Nucleotide::G, Nucleotide::A];

        let (ok, mismatches) = compare_with_tolerance(&a, &b, 1, 2);
        assert!(ok);
        assert_eq!(
            mismatches,
            Some(vec![(0, Nucleotide::T), (3, Nucleotide::A)])
        );
    }

    #[test]
    fn merge_left() {
        let seq_a = "ATTTCACATCTGTAATACTCTGTCTCCTTGTTATAATTTCATTTACTAGTTATAATTTATAATGCAAACTGGATTGCAGCCCCAGTGCCAGGACTCAAATTATCCCAGAAATATAGGAAAAAGATCAACTCACGGGGCTCCACGAAGAGTT";
        let seq_b = "GGCCTATTTCACATCTGTAATACTCTGTCTCCTTGTTATAATTTCATTTACTAGTTATAATTTATAATGCAAACTGGATTGCAGCCCCAGTGCCAGGACTCAAATTATCCCAGAAATATAGGAAAAAGATCAACTCACGGGGCTCCACGAA";

        let mut sequence_a = DNAString::new(seq_a.as_bytes().to_vec());
        let sequence_b = DNAString::new(seq_b.as_bytes().to_vec());

        match sequence_a.merge(&sequence_b, 20, 2, 1) {
            Ok(MergeResult::Left(_)) => { /* expected */ }
            Ok(other) => panic!("Unexpected merge result: {:?}", other),
            Err(err) => panic!("Merge failed: {:?}", err),
        }

        // Check that GGCCT was prepended
        let expected_prefix = DNAString::new("GGCCT".as_bytes().to_vec()).0;
        let actual_prefix = sequence_a.0.iter().take(5).cloned().collect::<Vec<_>>();
        assert_eq!(expected_prefix, actual_prefix);
    }

    #[test]
    fn merge_right() {
        let seq_a = "AATTGGCCAACCTTGGAAAG";
        let seq_b = "AAAGCTTTCGAC";

        let mut sequence_a = DNAString::new(seq_a.as_bytes().to_vec());
        let sequence_b = DNAString::new(seq_b.as_bytes().to_vec());

        match sequence_a.merge(&sequence_b, 4, 1, 1) {
            Ok(MergeResult::Right(_)) => { /* expected */ }
            Ok(other) => panic!("Unexpected merge result: {:?}", other),
            Err(err) => panic!("Merge failed: {:?}", err),
        }

        assert!(
            sequence_a.as_string().starts_with("AATTGGCC"),
            "Unexpected merge result"
        );
        assert!(
            sequence_a.as_string().ends_with("CTTTCGAC"),
            "Right extension failed"
        );
    }

    #[test]
    fn merge_inside() {
        let seq_a = "AAACCCGGGTTT";
        let seq_b = "CCCGGG";

        let mut sequence_a = DNAString::new(seq_a.as_bytes().to_vec());
        let sequence_b = DNAString::new(seq_b.as_bytes().to_vec());

        match sequence_a.merge(&sequence_b, 3, 0, 0) {
            Ok(MergeResult::Inside(_)) => { /* expected */ }
            Ok(other) => panic!("Unexpected merge result: {:?}", other),
            Err(err) => panic!("Merge failed: {:?}", err),
        }

        assert_eq!(sequence_a.as_string(), seq_a);
    }

    #[test]
    fn merge_outside() {
        let full_seq = "AAACCCGGGTTT";
        let sub_seq = "CCCGGG";

        let mut sequence_a = DNAString::new(sub_seq.as_bytes().to_vec());
        let sequence_b = DNAString::new(full_seq.as_bytes().to_vec());

        match sequence_a.merge(&sequence_b, 3, 0, 0) {
            Ok(MergeResult::Outside(_)) => { /* success */ }
            other => panic!("Unexpected result: {:?}", other),
        }

        assert_eq!(sequence_a.as_string(), full_seq);
    }

    #[test]
    fn merge_failure() {
        let seq_a = "AAAAAAA";
        let seq_b = "TTTTTTT";

        let mut sequence_a = DNAString::new(seq_a.as_bytes().to_vec());
        let sequence_b = DNAString::new(seq_b.as_bytes().to_vec());

        let result = sequence_a.merge(&sequence_b, 3, 1, 1);
        assert!(matches!(result, Err(MergeErrors::NoMerge)));
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
