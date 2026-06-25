// ==========================================
// Response to the Reviewers' Reports
// Manuscript: "Computational Études: A Spectral Approach"
// Addressee: Donna Chernyk, Senior Editor, Mathematics Editorial, Springer Nature
// Author: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
//         and Technology, Abu Dhabi, UAE)
// Created: June 25, 2026
// ==========================================

// ==========================================
// CONFIGURATION & TEMPLATE
// ==========================================

// Personal details
#let my-name = "Dr. Denys Dutykh"
#let my-title = "Acting Associate Dean of Graduate Studies"
#let my-college = "College of Computing and Mathematical Sciences"
#let my-email = "denys.dutykh@ku.ac.ae"
#let my-url = "https://www.denys-dutykh.com/"
#let my-uni = "Khalifa University of Science and Technology"
#let my-dpt = "Mathematics Department"
#let ku-blue = rgb("#005696")

// Manuscript title, reused throughout
#let book-title = "Computational Études: A Spectral Approach"

// Global page setup. The full letterhead is printed in the body of page 1 only;
// continuation pages carry a slim running header and a page-number footer.
#set page(
  paper: "a4",
  margin: (top: 3cm, bottom: 2.2cm, left: 2.5cm, right: 2.5cm),
  header: context {
    if counter(page).get().first() > 1 {
      set text(size: 8pt, fill: gray)
      [Response to the Referees · #emph(book-title)]
      v(-0.2em)
      line(length: 100%, stroke: 0.3pt + luma(200))
    }
  },
  footer: context {
    set text(size: 8pt, fill: gray)
    line(length: 100%, stroke: 0.3pt + luma(200))
    v(0.2em)
    grid(
      columns: (1fr, auto),
      align(left)[Khalifa University of Science and Technology · Abu Dhabi, UAE],
      align(right)[Page #counter(page).display("1") of #counter(page).final().first()],
    )
  },
)

// Global text and paragraph setup
#set text(font: "Linux Libertine O", size: 11pt, lang: "en")
#set par(justify: true, leading: 0.7em, spacing: 1.2em, first-line-indent: 0pt)

// Heading styling (used only in the enclosure)
#show heading.where(level: 1): it => block(below: 0.4em)[
  #set text(fill: ku-blue, size: 1.45em, weight: 700)
  #it.body
]
#show heading.where(level: 2): it => block(above: 1.4em, below: 0.6em)[
  #set text(fill: ku-blue, size: 1.12em, weight: 700)
  #it.body
]

// Helper: neutral grey-italic text marking an unfilled slot
#let pending(body) = text(fill: luma(130), style: "italic")[#body]

// Helper: the shaded box that holds a reviewer's verbatim report
#let referee-report(body) = block(
  width: 100%,
  fill: luma(246),
  inset: 11pt,
  stroke: (left: 2pt + ku-blue, rest: 0.5pt + luma(214)),
)[#body]

// Helper: one review section (report box followed by the author's reply).
// To complete the document, replace the two stub arguments of each
// call below with the verbatim report and the author's reply, respectively.
#let review-block(n, report, reply) = {
  heading(level: 2)[Review #n]
  // `sticky: true` keeps each bold label attached to the material that follows,
  // so a label is never stranded at the foot of a page once the real (longer)
  // reports and replies are pasted in.
  block(sticky: true)[#text(weight: "bold")[Referee's report.]]
  v(0.3em)
  referee-report(report)
  v(0.7em)
  block(sticky: true)[#text(weight: "bold")[Author's response.]]
  v(0.3em)
  reply
  v(1.4em)
}

// ==========================================
// PAGE 1: COVER LETTER TO THE EDITOR
// ==========================================

// Letterhead (first page only)
#grid(
  columns: (1fr, 1fr),
  align: (left + horizon, right + horizon),
  image("assets/KUlogo.png", width: 4.5cm),
  image("assets/CCMSLogo.png", width: 2.5cm),
)
#v(0.4em)
#line(length: 100%, stroke: 0.5pt + ku-blue)

#v(1.2em)

#align(right)[#datetime.today().display("[month repr:long] [day], [year]")]

#v(0.6em)

Donna Chernyk \
Senior Editor, Mathematics Editorial \
Springer Nature \
1650 Arch Street, Suite 2050 \
Philadelphia, PA 19103, USA

#v(1.1em)

#text(weight: "bold")[Subject: Response to the referees' reports on "#book-title"]

#v(0.7em)

Dear Ms. Chernyk,

Thank you for overseeing the review of my manuscript, #emph(book-title), and for keeping me so well informed about its progress. I am grateful to the referees for the time and care they have devoted to the text.

Please find enclosed my responses to the referees' reports, arranged in the order in which they were received and labelled Review 1 through Review 4. For each report I reproduce the Referee's comments in full and follow them with a single, considered response.

I have sought to engage each report as fully and as fairly as I can, revising the manuscript where the comments called for it. Should the referees or the Applied Mathematical Sciences series editors wish for any further clarification, I would be very glad to provide it, and I look forward to the next stage of the process.

#v(1.4em)

Sincerely,

#v(0.6em)

#grid(
  columns: (1fr, 1fr),
  gutter: 1em,
  align(left + horizon)[
    #text(weight: "bold")[#my-name] \
    #my-title \
    #my-dpt \
    #my-college \
    #my-uni \
    #link("mailto:" + my-email)[#text(fill: blue)[#my-email]] \
    #link(my-url)[#text(fill: blue)[#my-url]]
  ],
  align(right + horizon)[
    #image("assets/SignDD.png", width: 4cm)
  ],
)

#v(1.0em)

#text(style: "italic")[Enclosure: Responses to the four referees' reports.]

// ==========================================
// ENCLOSURE: RESPONSE TO THE REFEREES
// ==========================================

#pagebreak()

= Response to the Referees

#text(fill: luma(110), size: 10pt)[#emph(book-title) · #my-name]

#v(0.8em)

Throughout this enclosure, each Referee's report is set in the shaded box, with the author's response immediately beneath it. The reports are reproduced in full, save for the two most detailed ones, whose main points are summarised. References to pages, sections, and equations are to the revised manuscript.

#v(1.0em)

#review-block(1,
  [
    #emph[The first Referee provided an extensive and expert report, very positive in its overall judgement, with comments offered (in the Referee's words) "not critical" but for the author's consideration. Its main points are summarised here.]

    The Referee judges the book "very well written, excellently structured especially for (intended) educational purposes", and welcomes the inclusion of code in all three relevant languages (Python, MATLAB, and Julia), the exercises, and the emphasis on insight and illustration over formality and tedious error analysis, anticipating that it "will most certainly be appreciated by students". The substance of the report concerns the book's treatment of pseudospectral (PS) methods relative to finite differences (FD):

    - #strong[Central critique.] The manuscript compares PS chiefly against low-order (second-order) finite differences, whereas these are two ends of a continuum running through fourth and sixth order and beyond. The real disadvantages of PS, its limited geometric flexibility, the severe Chebyshev node clustering at the endpoints, the inability to refine locally, and the non-local character of spectral derivatives, are acknowledged but somewhat glossed over, while the recent literature that addresses them is largely absent from the discussion sections.

    - #strong[Terminology.] "Spectral" properly denotes the manipulation of expansion coefficients (the Galerkin approach), and "pseudospectral" the collocation methods that are the book's true subject; the distinction could be drawn earlier.

    - #strong[High-order finite differences.] Fourth-order FD is routine, and sixth order and above is of increasing importance (the subject of Fornberg's 2025 monograph); claims of PS superiority should be weighed against these. The label "myopic" for local stencils is questionable, since locality can be desirable, and the stated error gain of "a factor of four or eight, but rarely more" understates fourth-order behaviour.

    - #strong[Dimensionality and meshfree methods.] The deliberate one-dimensional focus sits awkwardly with the three-dimensional reality of applications, and the "A Broader Perspective" (§3.5) and "Modern Frontiers" (§16.29) sections omit radial basis function (RBF) and RBF-FD methods, which deliver spectral accuracy on scattered nodes in any dimension and recover sparse algebra. Several specific claims, that endpoint clustering is "essential", that node choice is "decisive for spectral accuracy", and that matrix-vector differentiation is special to PS, are not universally valid.

    - #strong[Further remarks.] The spectral advantage over high-order FD largely disappears for finitely-smooth data; dispersive errors and PS accuracy for non-smooth functions could be noted; boundary bordering may be simpler than deleting rows and columns; Gregory's quadrature deserves mention beside Newton--Cotes; and the "trapezoidal rule plus analytic correction" idea recurs in contour integration and in fractional-derivative computation. Back-references showing where each work is cited would help readers.
  ],
  [
    I am very grateful to the first Referee for so generous and so expert an assessment, and for comments offered, as the Referee graciously puts it, not in criticism but for the author's consideration. The report is plainly the work of a leading authority on high-order and meshfree methods, and engaging with it has improved the book considerably.

    The Referee's principal point is well made: the manuscript too often measured pseudospectral methods against low-order finite differences alone, and gave too little space to the modern high-order finite-difference and radial-basis-function literature that addresses the genuine limitations of the spectral approach. I agree without reservation. As scientists rather than advocates, we owe our students intellectual honesty above all; and while a foundational textbook cannot resolve what remain active research questions, it can present these trade-offs openly instead of glossing over them.

    The revised manuscript takes up the Referee's suggestions throughout. A terminological note now distinguishes, early in the Introduction, "spectral" methods in the modal or Galerkin sense from the "pseudospectral" collocation methods that are the book's true subject. The comparisons with finite differences have been rebalanced: the description of local stencils as "myopic" has been removed, with the explicit acknowledgement that locality is often a desirable property, and the effectiveness of high-order finite differences, sixth order and above, especially for the finitely-smooth data of real applications where Chebyshev clustering can be wasteful, is now stated plainly. The "A Broader Perspective" discussion and the "Modern Frontiers" material have been substantially expanded to treat radial basis function and RBF-FD methods, and the geometric flexibility they bring in two and three dimensions where tensor-product pseudospectral grids are weakest; the closing claim on the role of node choice has been qualified in the same spirit. Gregory's quadrature is now set beside Newton--Cotes, Fornberg's 2025 monograph has been added to the bibliography and cited at each of these points, and the Referee's further specific remarks have been worked into the relevant discussion sections.

    On the suggestion of back-references indicating where each work is cited: this is, regrettably, not yet reliably supported by the typesetting system used for the book, so I have been unable to include it in the present edition, though I would gladly add it once the feature matures. I thank the Referee once again for a report whose insight has measurably strengthened the manuscript, and for the warmth with which it received the book's pedagogical aims.
  ],
)

#review-block(2,
  [
    #emph[The Referee provided a detailed six-page report. It is warmly positive and recommends publication; its substantive content is a catalogue of minor corrections, whose main points are summarised here.]

    In overall terms the Referee judges the book "an ambitious, well-crafted graduate textbook on spectral methods", mathematically sound across the vast majority of its 535 pages, and singles out the three-language implementations, the historically deep literature surveys, and the étude format for particular praise. The issues raised are, with a single exception, typographical or structural:

    - #strong[One mathematical point (§7.5.2).] The Bernstein-ellipse parameter for the Witch of Agnesi is given as 1.5 rather than the golden ratio, with a consequent slip in the dependent claim about the rate of convergence.

    - #strong[Internal inconsistencies.] §1.6 describes a "dual-language" workflow and omits Julia, although the book is consistently tri-lingual; the "Literature Overview" headings alternate between sentence and title case across the chapters; and §4.11 and §8.8.1 give conflicting attributions for the discovery of the Padua points.

    - #strong[Table-of-contents gaps.] Computational Études 5.2 and 10.2 do not appear in the table of contents, and the §2.9 Exercises entry is missing.

    - #strong[Editorial and typographical items.] Two occurrences of "etude" lack the accent; an editorial self-correction was left in the printed text of Exercise 16.4; a small rounding discrepancy appears between the body and Table 2 in §3.3.2; and the Figure 3 caption states a period valid only for the demonstration parameters.

    - #strong[Smaller suggestions (non-errors).] A note on the editorial "we" in a single-authored book; a numbered equation and cross-reference for the Bernstein-ellipse formula; and an English gloss for the "Mise en Bouche" chapter title.

    The Referee concludes that one further proof-reading pass would bring the manuscript to publication-ready standard, and recommends the book warmly.
  ],
  [
    I thank the Referee most sincerely for an exceptionally careful and detailed reading of the manuscript, and for the generous overall assessment. The report is a model of constructive reviewing: each issue is pinpointed by section, quoted exactly, and paired with a suggested remedy, and I am grateful for the considerable time this must have taken.

    I am glad to report that every shortcoming the Referee identified has been corrected in the revised manuscript. The single mathematical point, the Bernstein-ellipse parameter for the Witch of Agnesi in §7.5.2, has been set to its correct value, the golden ratio $(1 + sqrt(5)) / 2 approx 1.618$, and the dependent convergence estimate adjusted to match. The remaining items, all of them typographical or structural, have likewise been resolved: the internal inconsistencies, the gaps in the table of contents, the missing accents, the stray editorial note in the exercises, and the minor numerical and caption discrepancies are all corrected, and the smaller stylistic suggestions have been taken up where they strengthen the book. I have also carried out, across the whole manuscript, the further proof-reading pass the Referee recommends. I am grateful for a review that has measurably improved the book, and for the warm recommendation that accompanies it.
  ],
)

#review-block(3,
  [
    I read the book from cover to cover and enjoyed it. There are two types of books on spectral methods: theory and practical. This book belongs to the second category, covering both elementary and more advanced topics. It should be understandable by students and scientists/engineers who do numerical computations.

    Here are my responses to your questions:

    #emph[How original is this work relative to existing books or papers in numerical analysis and computational PDEs?] \
    Readers can easily experiment with python or Matlab codes provided by the author. This makes it stand out among numerical PDE books.

    #emph[Are the computational experiments well designed, meaningful, and appropriately interpreted?] \
    Yes.

    #emph[Is the "computational études" format effective?] \
    Yes.

    #emph[Are there sections where clarity, motivation, or exposition could be significantly improved?] \
    No.

    #emph[Based on your overall assessment, would you recommend that Springer publish this manuscript?] \
    Yes.
  ],
  [
    I am very grateful to the Referee for reading the manuscript from cover to cover and for such an encouraging assessment. I am particularly glad that the Referee places the book in the practical rather than the purely theoretical category, and singles out the runnable codes that accompany every example as what sets it apart, since that hands-on, classroom-oriented character is exactly what the book sets out to achieve. I thank the Referee most warmly for the care taken over this report and for the recommendation that Springer publish the manuscript.
  ],
)

#review-block(4,
  [
    I had a quick read of the book draft. My overall impression is as follows:

    (i) The draft covers, or at least touches on, a wide range of topics, but the discussion is often too general and lacks depth and insight. Even for readers from the engineering community, a more focused and substantive treatment would be expected.

    (ii) The style and baseline of the book appear to be modeled on the books by Trefethen and J. P. Boyd. However, in my view, the current draft does not reach even half of the standard set by those works.

    (iii) The book offers limited value to researchers who intend to work on, or are already practicing, spectral methods.

    (iv) The draft contains many unnecessary descriptions and materials, which are distracting and weaken the focus of the presentation.

    Overall, my evaluation is that the book does not meet the standard expected for a Springer series.
  ],
  [
    I am grateful to the Referee for reading the draft and for stating their impression so plainly. Taken on its own terms, the assessment is a fair one, and the point of genuine difference is simply the readership the book is written for. Nearly every feature the Referee identifies as a shortcoming is a deliberate design choice, and, as it happens, one the manuscript already declares in its own front matter.

    The Referee's central observation, that the book "offers limited value to researchers who intend to work on, or are already practicing, spectral methods", is entirely correct, and it is so by design. #emph(book-title) is a graduate classroom textbook for interdisciplinary students meeting these methods for the first time; it is not a research monograph, nor a reference for specialists. This is not a distinction the reader is left to infer. The Preface carries a section headed "What This Book Is, and What It Is Not", and the Introduction opens with the same declaration, stating plainly that the book is "written to be read, and run, from cover to cover by a student who does not yet have a spectral-methods specialist looking over their shoulder", and pointing the specialist instead to the standard references for the complete theory and the current state of the art. Read with that audience in mind, the remaining points largely answer themselves.

    The Referee is right, too, that the book takes its baseline from Trefethen, Boyd, and Fornberg. This is an openly acknowledged homage, set out in the Preface, and I hold those works to be masterpieces that the present text does not pretend to rival in depth. What has changed since they were written is the setting in which a student learns. Their delivery rests on an older computational stack, and when an instructor today places Fortran code in the hands of an interdisciplinary doctoral student, the pedagogical transmission stops almost at once. The aim of this book is to carry the same philosophy natively into the tools today's students use, with every example worked in Python, Julia, and MATLAB, so that the reader runs and adapts the code rather than transcribes it.

    The breadth that the Referee reads as a want of depth is, in the same way, the organising principle of the book, which favours horizontal reach over vertical exhaustion. The material was developed and tested in live graduate classrooms with mixed cohorts, and its purpose is to give a student in fluid dynamics, a student in solid mechanics, and a student in electromagnetics a common and workable foundation, from which each can take the one technique their own thesis requires. A treatment narrowed and deepened for the specialist would serve that shared purpose less well, not better.

    The same consideration answers the concern about unnecessary or distracting material. Much of what reads as redundant from the vantage point of expertise is precisely the explanation that allows a first-year graduate student, struggling to make a differentiation matrix converge, to succeed without a specialist at their side. The book dwells on these small, practical matters deliberately, and for its intended reader they are often the most valuable pages in the chapter.

    I therefore stand by the architecture and the intent of the manuscript. That it could be read as a research monograph that falls short, rather than as the classroom text it declares itself to be, tells me its statement of purpose must be impossible to overlook; it now stands plainly at the front of both the Preface and the Introduction, so that no future reader need make the same, entirely understandable, mistake. Assessed as what it is, a modern and broadly adoptable graduate textbook that a student can read from cover to cover and come away able to use, I believe it meets the standard that purpose sets, and I hope it may be judged in that light.
  ],
)
