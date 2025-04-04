# merge_sequences
![example workflow](https://github.com/Dr-TSteimle/merge_sequences/actions/workflows/rust.yml/badge.svg)

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
merge_sequences  = { git = "https://github.com/Dr-TSteimle/merge_sequences.git" }
```
