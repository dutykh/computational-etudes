// style/template.typ

// Import the droplet package for proper drop caps
#import "@preview/droplet:0.3.1": dropcap as droplet-dropcap
// Import indenta for automatic first-paragraph indent handling
#import "@preview/indenta:0.0.3": fix-indent
// Import codly for beautiful code blocks
#import "@preview/codly:1.3.0": *
#import "@preview/codly-languages:0.1.1": *

// --- DROP CAP FUNCTION ---
// Wrapper around droplet package with our styling
#let dropcap(body) = {
  droplet-dropcap(
    height: 3,
    gap: 4pt,
    overhang: 0pt,
    font: "New Computer Modern",
    weight: "bold",
    fill: rgb(20, 45, 110),
    body
  )
}


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
  set text(font: "New Computer Modern", lang: "en", size: 11pt)
  set par(
    justify: true,
    first-line-indent: 1.5em,
    leading: 0.65em,
    spacing: 0.55em,
  )

  // --- EM-DASH HANDLING ---
  show "---": [ #h(0.1em)—#h(0.1em) ]

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

  // --- HEADING STYLING ---
  set heading(numbering: "1.1")

  // Level 1 Heading (Chapter) design
  show heading.where(level: 1): it => {
    pagebreak(weak: true) // Always start chapters on a new page
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
        #counter(heading).display(it.numbering) #h(0.4em) #it.body
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

  // --- FRONT MATTER (Preface, TOC) ---
  set page(
    paper: "a4",
    margin: (inside: 3cm, outside: 2cm, y: 2.5cm),
    numbering: "i",
  )
  counter(page).update(1)

  // Table of Contents
  {
    set par(leading: 1.2em, first-line-indent: 0em) // Spacing between lines

    // Make TOC entries use the Navy color
    show outline.entry: it => text(fill: navy, it)

    // Make all hyperlinks use the Navy color
    show link: it => text(fill: navy, it)

    outline(depth: 2, indent: auto)
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

  // Apply fix-indent after all other show rules for automatic first-paragraph handling
  show: fix-indent()

  body
}
