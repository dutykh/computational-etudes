// textbook/chapters/appendix_e_hints.typ
// Appendix E: Hints and Solutions to Selected Exercises
// Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
// Last modified: June 2026

#import "../styles/template.typ": dropcap, chapter-abstract

= Hints and Solutions to Selected Exercises <appendix-hints>

#chapter-abstract(keywords: [Exercises · Hints · Solutions])[
This appendix gathers hints and short solutions for a selection of the harder exercises spread across the chapters, principally the project-style problems and the more demanding computational ones. Each entry is keyed to its exercise number and links back to the statement. The remaining exercises are left to the reader, in keeping with the active-learning spirit of the computational études.
]

#dropcap[Spectral methods reward the reader who computes alongside the text, and the exercises are written to be attempted before they are read about.] The notes below are deliberately brief. They name the key idea, point to the relevant equation or étude, or take the first step of a derivation, rather than working every detail; a reader who has already wrestled with a problem usually needs no more than a sentence to get unstuck. Where an exercise asks for a numerical experiment, the hint states the expected qualitative outcome and the diagnostic that confirms it.

#let navy = rgb(20, 45, 110)

// The hints are authored beside their exercises in the chapter sources via the
// deferred `#hint-for(<label>)[...]` / `#solution-for(<label>)[...]` markers, and
// collected here in document order, each keyed to its auto-numbered exercise.
#context {
  let entries = query(<solution-entry>)
  for s in entries {
    let v = s.value
    let lead = if v.kind == "hint" { [Hint] } else { [Solution] }
    block(spacing: 0.85em)[
      #text(weight: "semibold", fill: navy)[#lead to #ref(v.ex).] #h(0.3em) #v.body
    ]
  }
}
