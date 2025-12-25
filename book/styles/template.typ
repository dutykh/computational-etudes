// style/template.typ

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
  
  set page(
    paper: "a4",
    margin: (inside: 3cm, outside: 2cm, y: 2.5cm),
    numbering: "1",
  )

  // Use New Computer Modern for that "classic math" look
  set text(font: "New Computer Modern", lang: "en", size: 11pt)
  set par(justify: true)

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
        #text(size: 2.2em, weight: 700, fill: navy, it.body)
      ]
    )
    v(1.5cm) // Space after header
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
          set par(justify: false) // Prevents "Com-pu-ta-tion-al" spacing
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
  set page(numbering: "i")
  counter(page).update(1)
  
  // Table of Contents
  {
    set par(leading: 1.2em) // Spacing between lines
    
    // Make TOC entries use the Navy color
    show outline.entry: it => text(fill: navy, it)
    
    // Make all hyperlinks use the Navy color
    show link: it => text(fill: navy, it)
    
    outline(depth: 2, indent: auto)
  }

  // --- MAIN CONTENT START ---
  set page(numbering: "1")
  counter(page).update(1)

  body
}
