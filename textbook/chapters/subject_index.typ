// textbook/chapters/subject_index.typ
// Subject Index — collected automatically from #idx("term") markers.
// Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

#heading(level: 1, numbering: none)[Subject Index]

#let _navy = rgb(20, 45, 110)

The page references below collect every `#idx("term")` marker planted in the body of the volume.  Pages are absolute (in the book's main numbering, not the per-chapter numbering).  Multiple page references for the same term are listed in ascending order, separated by commas.  This index is *selective rather than exhaustive*: terms that appear only in passing are left to the reader; terms that name a recurring tool, method, equation, or named theorem are indexed at *every* page on which they substantively appear.

#v(0.6em)

#context {
  // Collect every <index-entry> marker.
  let entries = query(<index-entry>)

  // Map: term -> sorted list of unique page numbers.
  let terms = (:)
  for e in entries {
    let term = e.value
    let p = counter(page).at(e.location()).first()
    if term in terms {
      if p not in terms.at(term) {
        terms.insert(term, terms.at(term) + (p,))
      }
    } else {
      terms.insert(term, (p,))
    }
  }

  // Sort the term keys alphabetically (case-insensitive).
  let sorted-terms = terms.keys().sorted(key: t => lower(t))

  if sorted-terms.len() == 0 {
    [_The index is empty: no `#idx(...)` markers have been planted yet._]
  } else {
    // Group terms by leading letter for visual grouping.
    let groups = (:)
    for t in sorted-terms {
      let first = upper(t.at(0))
      // Treat anything not starting with A--Z as "Other".
      if not ("A" <= first and first <= "Z") {
        first = "Other"
      }
      if first in groups {
        groups.insert(first, groups.at(first) + (t,))
      } else {
        groups.insert(first, (t,))
      }
    }

    // Render the two-column index.
    set par(justify: false, first-line-indent: 0em, leading: 0.55em)
    show: columns.with(2, gutter: 1.2em)

    for letter in groups.keys().sorted() {
      block(below: 0.4em)[
        #text(size: 1.1em, weight: "bold", fill: _navy)[#letter]
      ]
      for t in groups.at(letter).sorted() {
        let pages = terms.at(t).sorted()
        block(below: 0.25em, breakable: false)[
          #t, #pages.map(p => str(p)).join(", ")
        ]
      }
      v(0.5em)
    }
  }
}
