// style/template.typ

// Import the droplet package for proper drop caps
#import "@preview/droplet:0.3.1": dropcap as droplet-dropcap
// First-paragraph indentation uses Typst's native `first-line-indent` with
// `all: true` (set in `project` below): every main-flow paragraph is indented,
// while block interiors (the drop cap, figure captions, the left-rule boxes)
// stay flush. This works uniformly across the `#include`d chapters, which a
// top-level show rule cannot reach.
// Import codly for beautiful code blocks
#import "@preview/codly:1.3.0": *
#import "@preview/codly-languages:0.1.1": *
// Import zero for scientific number formatting
#import "@preview/zero:0.6.1": num, format-table, set-num, set-group

// --- DROP CAP FUNCTION ---
// Wrapper around droplet package with our styling
#let dropcap(body) = {
  // Wrapped in a block() so the `fix-indent` show rule classifies the drop cap
  // as a block, and therefore indents the paragraph that follows it.
  block(breakable: true, above: 0em, below: 0.55em)[
    #droplet-dropcap(
      height: 3,
      gap: 4pt,
      overhang: 0pt,
      font: ("Linux Libertine O"),
      weight: "bold",
      fill: rgb(20, 45, 110),
      body,
    )
  ]
}

// --- MATHEMATICAL CONSTANTS ---
// Upright imaginary unit (like \mathrm{i} in LaTeX)
// Defined as a standard upright math atom, just like "dif" for proper spacing
#let ii = math.upright("i")

// --- SUBJECT-INDEX MARKER ---
// `#idx("term")` plants an invisible labelled `metadata` element at the
// cursor.  The Subject Index page in the back matter queries every
// `<index-entry>` label, sorts by term, groups duplicates with their
// page numbers, and emits a two-column index in the same navy-rule
// style as the notation glossary.
#let idx(term) = [#metadata(term)<index-entry>]

// --- COMPUTATIONAL-ÉTUDE CONCLUSION BLOCK ---
// Closing paragraph for a Computational Étude.  Renders a thin sky-blue
// left rule with a bold-navy `Takeaway.` lead-in.  Lighter than the
// navy-bordered Principle / Definition / Rule-of-Thumb boxes so that
// conclusions remain subordinate to named mathematical results.
#let etude-conclusion(body) = {
  let navy = rgb(20, 45, 110)
  let sky = rgb(120, 150, 210)
  v(0.6em)
  block(
    stroke: (left: 1.2pt + sky),
    inset: (left: 12pt, top: 4pt, bottom: 4pt, right: 0pt),
    spacing: 0.65em,
  )[
    #text(weight: "semibold", fill: navy)[Takeaway.] #h(0.4em) #body
  ]
  v(0.4em)
}

// --- CHAPTER ABSTRACT BLOCK ---
// Springer-style chapter abstract (+ optional keywords) for SpringerLink.
// Rendered as a thin sky-blue left rule with a bold-navy "Abstract" label on
// its own line, echoing the lighter `etude-conclusion` style so the abstract
// reads as introductory matter rather than a named mathematical result.
// Placed immediately after a chapter's level-1 heading and before the opening
// drop cap.  `keywords` is optional inline content (e.g. a `·`-separated run).
#let chapter-abstract(body, keywords: none) = {
  let navy = rgb(20, 45, 110)
  let sky = rgb(120, 150, 210)
  block(
    width: 100%,
    stroke: (left: 1.2pt + sky),
    inset: (left: 12pt, top: 4pt, bottom: 6pt, right: 0pt),
    spacing: 0.8em,
  )[
    #set par(first-line-indent: 0em, justify: true, leading: 0.6em)
    #set text(size: 0.95em)
    #text(weight: "semibold", fill: navy)[Abstract]
    #v(0.45em)
    #body
    #if keywords != none {
      v(0.5em)
      [#text(weight: "semibold", fill: navy)[Keywords] #h(0.5em) #keywords]
    }
  ]
  v(0.8em)
}

// --- EXERCISE ENVIRONMENT ---
// A numbered, cross-referenceable exercise in the house left-rule style.
// Implemented as a `figure` with a custom `kind: "exercise"` so that
// per-chapter numbering is INHERITED from the document-wide
// `set figure(numbering: ...)` (body chapters get 8.1, appendices get A.1 with
// no extra code), `@ex-...` cross-references resolve to "Exercise 8.1" for
// free, and the set can be collected into a "List of Exercises" exactly like
// the Computational-Étude list.  The optional short `title` rides in the figure
// `caption`; the show rule reads it back as `it.caption.body` and it is never
// realised as a caption.
#let exercise(title: none, body) = figure(
  kind: "exercise",
  supplement: [Exercise],
  caption: title,
  body,
)

// --- INLINE HINT / SOLUTION BLOCKS ---
// Lighter than the named-result boxes, echoing `etude-conclusion`.  These
// render at the call site; use the deferred forms below to gather answers in a
// back-matter appendix instead.
#let hint(body) = {
  let teal = rgb(20, 130, 130)
  v(0.4em)
  block(
    stroke: (left: 1.0pt + teal),
    inset: (left: 12pt, top: 3pt, bottom: 3pt, right: 0pt),
    spacing: 0.6em,
  )[#text(weight: "semibold", fill: teal)[Hint.] #h(0.4em) #body]
  v(0.3em)
}

#let solution(body) = {
  let navy = rgb(20, 45, 110)
  let sky = rgb(120, 150, 210)
  v(0.4em)
  block(
    stroke: (left: 1.0pt + sky),
    inset: (left: 12pt, top: 3pt, bottom: 3pt, right: 0pt),
    spacing: 0.6em,
  )[#text(weight: "semibold", fill: navy)[Solution.] #h(0.4em) #body]
  v(0.3em)
}

// --- DEFERRED HINTS / SOLUTIONS ---
// `#hint-for(<ex-bvp-greens>)[ ... ]` attaches a hint to a labelled exercise.
// Nothing renders at the call site; the "Hints and Solutions" appendix queries
// every <solution-entry> and prints it keyed to its exercise number via `ref`.
#let solution-for(ex, body) = [#metadata((kind: "solution", ex: ex, body: body))<solution-entry>]
#let hint-for(ex, body) = [#metadata((kind: "hint", ex: ex, body: body))<solution-entry>]

#let project(
  title: "",
  subtitle: "",
  author: "",
  affiliation: [],
  date: none,
  body,
) = {
  // --- COLORS ---
  let navy = rgb(20, 45, 110)
  let sky = rgb(120, 150, 210)

  // --- SCIENTIFIC NUMBER FORMATTING (zero package) ---
  set-num(product: sym.times)  // Use × for scientific notation
  set-group(threshold: 5)       // Only group if >= 5 digits

  // --- GLOBAL DOCUMENT SETTINGS ---
  set document(author: author, title: title)

  show link: set text(fill: navy)

  // --- CODE BLOCK STYLING (codly) ---
  show: codly-init.with()
  codly(
    languages: codly-languages,
    zebra-fill: none,
    stroke: 0.5pt + luma(200),
    fill: luma(250),
    radius: 3pt,
    number-format: (n) => text(fill: luma(140), size: 0.85em, str(n)),
  )

  // --- PARAGRAPH TYPOGRAPHY ---
  set text(font: ("Linux Libertine O", "New Computer Modern"), lang: "en", size: 11pt)
  set par(
    justify: true,
    first-line-indent: (amount: 1.5em, all: true),
    leading: 0.65em,
    spacing: 0.55em,
  )

  // --- LIST STYLING ---
  set list(
    indent: 1.5em,
    body-indent: 0.5em,
    spacing: 0.65em,
  )

  set enum(
    indent: 1.5em,
    body-indent: 0.5em,
    spacing: 0.65em,
  )

  // Add vertical space around lists
  show list: it => {
    v(0.6em)
    it
    v(0.6em)
  }

  show enum: it => {
    v(0.6em)
    it
    v(0.6em)
  }

  // --- FIGURE AND TABLE STYLING ---
  // Add breathing room around figures and tables
  show figure: it => {
    v(1.5em)
    it
    v(1.5em)
  }

  // Style captions: smaller, italic body, bold navy prefix
  show figure.caption: it => {
    set text(size: 0.92em)
    set par(justify: true, leading: 0.5em)
    v(0.6em)
    block(width: 90%, inset: (x: 5%))[
      #text(weight: "semibold", fill: navy)[#it.supplement #context it.counter.display(it.numbering):]
      #h(0.3em)
      #emph[#it.body]
    ]
  }

  // --- EXERCISE BLOCK STYLING ---
  // Mirrors the `etude-conclusion` left-rule house style.  Being more specific
  // than, and defined after, the generic `show figure` rule above, this rule
  // fully replaces default figure rendering for exercises and supplies its own
  // vertical spacing.  `breakable: true` lets a long multi-part exercise split
  // across a page boundary.
  show figure.where(kind: "exercise"): it => {
    v(0.7em)
    block(
      width: 100%,
      breakable: true,
      stroke: (left: 1.2pt + sky),
      inset: (left: 12pt, top: 5pt, bottom: 5pt, right: 0pt),
      spacing: 0.65em,
    )[
      #text(weight: "semibold", fill: navy)[Exercise #context it.counter.display(it.numbering)]#if it.caption != none [ (#emph(it.caption.body)).] else [.] #h(0.4em) #it.body
    ]
    v(0.5em)
  }

  // --- BIBLIOGRAPHY STYLING ---
  // Use a real heading so it appears in the outline/TOC.
  // The level-1 heading show rule handles the visual styling (pagebreak, navy bar).
  // Setting title: none suppresses the built-in bibliography heading.
  show bibliography: it => {
    heading(level: 1, numbering: none)[Bibliography]
    it
  }

  set bibliography(title: none, style: "numeric-alphabetical.csl")

  // --- HEADING STYLING ---
  set heading(numbering: "1.1")

  // --- PER-CHAPTER EQUATION NUMBERING ---
  // Equation labels are formatted as (chapter.local), e.g. (20.3).  The local
  // counter is reset at the start of every top-level heading.
  set math.equation(numbering: n => {
    let chapters = counter(heading).get()
    if chapters.len() > 0 and chapters.first() > 0 {
      numbering("(1.1)", chapters.first(), n)
    } else {
      numbering("(1)", n)
    }
  })

  // --- PER-CHAPTER FIGURE & TABLE NUMBERING ---
  // Figure and table labels are formatted as chapter.local, e.g. 20.3.  Both
  // kind-specific counters are reset at the start of every top-level heading.
  set figure(numbering: n => {
    let chapters = counter(heading).get()
    if chapters.len() > 0 and chapters.first() > 0 {
      numbering("1.1", chapters.first(), n)
    } else {
      numbering("1", n)
    }
  })
  // Level 1 Heading (Chapter) design.  We ALSO reset the math equation, figure,
  // and table counters here so that per-chapter labels (e.g. (20.3), Fig 20.1)
  // start fresh on every new chapter.
  show heading.where(level: 1): it => {
    pagebreak(weak: true) // Always start chapters on a new page
    counter(math.equation).update(0)
    counter(figure.where(kind: image)).update(0)
    counter(figure.where(kind: table)).update(0)
    counter(figure.where(kind: "exercise")).update(0)
    let number = if it.numbering == none { none } else { counter(heading).display(it.numbering) }

    v(3cm) // Vertical space at top of chapter

    // The Chapter visual container
    grid(
      columns: (6mm, 1fr),
      gutter: 5mm,
      // The blue sidebar for the header
      rect(width: 100%, height: 100%, fill: navy),
      block[
        #if number != none [
           #text(size: 1.2em, weight: "bold", fill: sky)[CHAPTER #number]
           #v(0.3em)
        ]
        #set text(hyphenate: false)
        #text(size: 2.0em, weight: 700, fill: navy, it.body)
      ]
    )
    v(1.5cm) // Space after header
  }

  // Level 2 Heading with decorative rule
  show heading.where(level: 2): it => {
    v(1.8em)
    block(below: 0.8em)[
      #text(size: 1.35em, weight: "semibold", fill: navy)[
        #if it.numbering != none [
          #counter(heading).display(it.numbering) #h(0.4em)
        ]
        #it.body
      ]
      #v(0.3em)
      #line(length: 2cm, stroke: 0.75pt + sky)
    ]
  }

  // Level 3 Heading (italic style)
  show heading.where(level: 3): it => {
    v(1.2em)
    block(below: 1.0em)[
      #text(size: 1.1em, weight: "medium", style: "italic", fill: navy)[
        #if it.numbering != none [
          #counter(heading).display(it.numbering) #h(0.3em)
        ]
        #it.body
      ]
    ]
  }

  // Level 4 Heading (bold, no numbering)
  show heading.where(level: 4): it => {
    v(1em)
    block(below: 0.8em)[
      #text(size: 1em, weight: "semibold", fill: navy)[#it.body]
    ]
  }

  // --- TITLE PAGE ---
  // We use a dedicated page block with specific margins
  page(margin: (top: 0cm, bottom: 0cm, left: 0cm, right: 0cm), numbering: none)[

    // 1. The Sidebar Strip (Placed absolutely to ensure full height)
    #place(left, rect(width: 3cm, height: 100%, fill: navy))

    // 2. The Content Container (Offset to right of the strip)
    #place(
      top + left,
      dx: 3cm, // Start after the blue strip
      dy: 0cm,
      block(
        width: 100% - 3cm, // Remaining width
        height: 100%,
        inset: (x: 1.5cm, y: 3cm), // Internal padding
        breakable: false,
        {
          // Title Page Typography Settings
          set par(justify: false, first-line-indent: 0em) // Prevents spacing issues
          set text(hyphenate: false) // Prevents word breaking in titles

          // Series Title
          align(left)[
            #text(0.9em, fill: gray, weight: 600, tracking: 0.1em)[
              NOTES AND STUDIES IN SPECTRAL METHODS
            ]
            #v(0.2em)
            #line(length: 4cm, stroke: 1pt + sky)
          ]

          v(2fr) // Flexible space

          // Main Title
          align(left)[
            // Fixed leading by setting it on the paragraph, not the text
            #set par(leading: 0.3em)
            #text(3.5em, weight: 800, fill: navy, title)
            #v(0.5em)
            #text(1.8em, weight: 500, style: "italic", fill: navy.lighten(20%), subtitle)
          ]

          v(2fr)

          // Author & Affiliation
          align(left)[
            #text(1.4em, weight: 700, fill: black, author)
            #v(1em)
            // Affiliation block with smaller, lighter text
            #block(width: 80%)[
              #set text(size: 0.9em, fill: luma(40%), style: "italic")
              #affiliation
            ]
          ]

          v(1fr)

          // Date
          align(bottom + left)[
             #text(1em, weight: 600, fill: navy)[
               #if date != none { date } else { datetime.today().display("[year]") }
             ]
          ]

          v(1cm)
        }
      )
    )
  ]

  // --- DEDICATION ---
  // A single, unnumbered page placed immediately after the title page,
  // before the front-matter Roman numbering begins.  Vertically centred,
  // italic, with generous whitespace.
  page(margin: (top: 0cm, bottom: 0cm, left: 0cm, right: 0cm), numbering: none)[
    #v(1fr)
    #align(center)[
      #set par(justify: false, first-line-indent: 0em, leading: 1em)
      #set text(size: 1.15em, style: "italic", fill: navy)
      To Katya, my wife,
      #v(0.6em)
      and to our sons Michel and Nicolas.
    ]
    #v(2fr)
  ]

  // --- FRONT MATTER (Preface, TOC) ---
  set page(
    paper: "a4",
    margin: (inside: 3cm, outside: 2cm, y: 2.5cm),
    numbering: "i",
  )
  counter(page).update(1)

  // Table of Contents
  {
    set par(leading: 0.65em, first-line-indent: 0em) // Spacing between lines

    // Make TOC entries use the Navy color
    show outline.entry: it => text(fill: navy, it)

    // Make all hyperlinks use the Navy color
    show link: it => text(fill: navy, it)

    outline(depth: 2, indent: auto)
  }

  // --- LIST OF COMPUTATIONAL ÉTUDES ---
  // The étude is the central pedagogical unit of the book.  We collect every
  // heading whose body text begins with "Computational Étude" (regardless of
  // whether it sits at level 2 or level 3 in the chapter file) and emit a
  // dotted-leader outline keyed on its absolute page number.
  pagebreak()
  {
    set par(leading: 0.65em, first-line-indent: 0em)
    show link: it => text(fill: navy, it)
    heading(level: 1, numbering: none)[List of Computational Études]

    context {
      // Filter on heading level 2 or 3 (études) AND body containing
      // "Computational Étude".  The level filter excludes the LoCe heading
      // itself (level 1), and the text filter selects only the étude
      // headings rather than every section in the chapter.
      let etudes = query(heading).filter(h =>
        h.level >= 2 and repr(h.body).contains("Computational Étude")
      )
      for h in etudes {
        let pn = counter(page).at(h.location()).first()
        link(h.location())[
          #text(fill: navy)[
            #h.body
            #box(width: 1fr, repeat[ . ])
            #pn
          ]
        ]
        linebreak()
      }
    }
  }

  // --- MAIN CONTENT START ---
  set page(
    paper: "a4",
    margin: (inside: 3cm, outside: 2cm, top: 2.5cm, bottom: 2.5cm),
    numbering: "1",
    header: context {
      // Get the current chapter
      let chapters = query(selector(heading.where(level: 1)).before(here()))
      if chapters.len() > 0 {
        let current-chapter = chapters.last()
        // Don't show header on chapter opening pages
        let chapter-pages = query(selector(heading.where(level: 1)))
        let on-chapter-page = chapter-pages.any(h => {
          let h-loc = h.location()
          let here-loc = here()
          h-loc.page() == here-loc.page()
        })

        if not on-chapter-page {
          set text(size: 9pt, fill: luma(120))
          grid(
            columns: (1fr, auto),
            align(left, smallcaps(current-chapter.body)),
            align(right, counter(page).display()),
          )
          v(-0.3em)
          line(length: 100%, stroke: 0.5pt + luma(200))
        }
      }
    },
    header-ascent: 40%,
  )
  counter(page).update(1)

  body
}
