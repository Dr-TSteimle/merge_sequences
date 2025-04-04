# merge_sequences

<<<<<<< HEAD
**Approximate merging of DNA sequences using k-mer overlaps and mismatch-tolerant alignment**  
_A Rust library by Thomas Steimlé_

---

## 🔬 What is it?

`merge_sequences` provides a fast, accurate, and configurable way to **merge partially overlapping DNA sequences**, accounting for errors or mismatches. It's designed for use cases like:

- Read clustering
- Contig extension
- De novo assembly graph simplification

Merging is based on **shared k-mer overlaps**, with support for mismatch tolerance and tracking of alignment regions.

---

## ⚙️ Features

- ✅ K-mer based overlap detection
- ✅ Mismatch-tolerant alignment (`max_mismatches`, `max_consecutive_mismatches`)
- ✅ Merge classification: `Left`, `Right`, `Inside`, `Outside`
- ✅ DNA abstraction (`DNAString`, `Nucleotide`)
- ✅ Detailed mismatch reporting

---

## 📦 Installation

Add the following to your `Cargo.toml`:

```toml
[dependencies]
merge_sequences = "0.1"

=======
![example workflow](https://github.com/Dr-TSteimle/merge_sequences/actions/workflows/rust.yml/badge.svg)
>>>>>>> d9116841f8df0fa3aea9ed726d232833ebb42857
