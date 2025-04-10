![example workflow](https://github.com/Dr-TSteimle/merge_sequences/actions/workflows/rust.yml/badge.svg)

# merge_sequences
<<<<<<< HEAD

**Approximate merging of DNA sequences using k-mer overlaps and mismatch-tolerant alignment**  

---
=======
![example workflow](https://github.com/Dr-TSteimle/merge_sequences/actions/workflows/rust.yml/badge.svg)
>>>>>>> a8db3c78e633ed654e6721dacf257f9831831368

## 🔬 What is it?

`merge_sequences` provides a fast, accurate, and configurable way to **merge partially overlapping DNA sequences**, accounting for errors or mismatches. It's designed for use cases like:

- Read clustering
- Contig extension

Merging is based on **shared k-mer overlaps**, with support for mismatch tolerance and tracking of alignment regions.

---

## ⚙️ Features

- ✅ K-mer based overlap detection
- ✅ Mismatch-tolerant alignment (`max_mismatches`, `max_consecutive_mismatches`)
- ✅ Merge classification: `Left`, `Right`, `Inside`, `Outside`

---

## 📦 Installation

Add the following to your `Cargo.toml`:

```toml
[dependencies]
<<<<<<< HEAD
merge_sequences = "0.1"
```

=======
merge_sequences  = { git = "https://github.com/Dr-TSteimle/merge_sequences.git" }
```
>>>>>>> a8db3c78e633ed654e6721dacf257f9831831368
