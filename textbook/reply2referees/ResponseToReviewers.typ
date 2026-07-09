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

Please find enclosed my responses to the referees' reports, arranged in the order in which they were received and labelled Review 1 through Review 6. For each report I reproduce the Referee's comments in full and follow them with a single, considered response.

I have sought to engage each report as fully and as fairly as I can, revising the manuscript where the comments called for it. Following your guidance, and that of the two most recent reports, I have also undertaken a substantial shortening of the manuscript, whose specific measures are set out in the responses that follow. Should the referees or the Texts in Applied Mathematics series editors wish for any further clarification, I would be very glad to provide it, and I look forward to the next stage of the process.

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

#text(style: "italic")[Enclosure: Responses to the six referees' reports.]

// ==========================================
// ENCLOSURE: RESPONSE TO THE REFEREES
// ==========================================

#pagebreak()

= Response to the Referees

#text(fill: luma(110), size: 10pt)[#emph(book-title) · #my-name]

#v(0.8em)

Throughout this enclosure, each Referee's report is set in the shaded box, with the author's response immediately beneath it. The reports are reproduced in full, save for the four most detailed ones, whose main points are summarised. References to pages, sections, and equations are to the revised manuscript.

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

#review-block(5,
  [
    #emph[The fifth Referee, who has taught spectral methods over many years, provided a four-page report recommending against publication in the present form. It is the most critical of the five, and its main points are summarised here under the Referee's own headings.]

    #strong[Overview and assessment.] The Referee welcomes the project's aim, a modern text that synthesises and updates the two classics of Trefethen (2000) and Boyd (2000), combining the accessibility of the former with the breadth and practicality of the latter, and calls the idea "a good one". The recommendation is nonetheless against publication in the current form, pending a careful reworking of the bibliography, a more balanced treatment of modern developments, and improvements to the exposition; the stated reservation is one of confidence in the manuscript's "scholarly judgment, organization or technical precision".

    #strong[Main concern: scholarship and bibliography.] The Referee's principal concern is that the engagement with advanced and modern topics does not always distinguish established developments from more marginal or preliminary ones. Several references are judged to carry few citations or to appear in obscure venues, and to be presented without sufficient context explaining why they were selected or how they relate to the established literature; the worry is both that undue prominence is given to low-impact work and that important established references may be underemphasised. Fourier continuation and FFT-based methods for non-periodic problems are cited as a substantial body of work that appears to be missing.

    #strong[Scope and audience.] At over 750 pages the manuscript is judged too long for a typical course, yet, in the Referee's view, not authoritative enough to serve as a reference for graduate students or researchers; the intended audience is therefore felt to be unclear.

    #strong[Presentation and exposition.] The Referee finds the early exposition, in places, choppy or imprecise: Chebyshev points are used and mentioned before they are formally defined, and the index points to a mention rather than to the definition; Legendre nodes likewise; the barycentric interpolation formula is treated very briefly beside the lengthy development of Lebesgue constants, though it is the more practically useful tool; the Bernstein-ellipse notation is used before it is defined; the geometric-convergence rate $O(rho^(-N))$ is stated without the algebraic factor that can accompany a singularity lying on the ellipse; and Computational Étude 4.1, which asks whether the endpoint clustering of Chebyshev points is essential, studies only uniformly random nodes. The Referee observes that the real difficulty is independent sampling, which places points arbitrarily close together, and notes that a determinantal point process with the same arcsine marginal but added repulsion yields far better Lebesgue constants.
  ],
  [
    I am grateful to the fifth Referee, whose long experience of teaching spectral methods is evident throughout, for a close and candid reading. This is the most critical of the five reports, and I have treated it accordingly. Where its specific criticisms are well founded I have adopted them, and the manuscript is the better for it; where they rest on a view of the book's purpose that the book itself does not claim, I hope to explain the design plainly. I take the points in the order the Referee raises them, beginning with those I have acted upon.

    Several of the Referee's observations on the early exposition were simply correct, and I have made the corresponding revisions. The chapter on the geometry of nodes now defines the Bernstein ellipse $cal(E)_rho$ explicitly, as the image of the circle $abs(z) = rho$ under the Joukowski map, before the notation is relied upon. The geometric-convergence estimate is now stated with the caveat the Referee rightly asks for: the clean $O(rho^(-N))$ bound holds for a function analytic strictly inside $cal(E)_rho$, whereas a singularity lying exactly on the ellipse can contribute a factor growing algebraically in $N$, as Boyd discusses. The barycentric interpolation formula has been given the fuller treatment its practical importance deserves, its derivation, its scale-invariance, and Higham's stability guarantee now spelled out, since the Referee is right that it is more useful to most readers than the neighbouring Lebesgue-constant bounds. And the ordering complaints have been addressed at their root: the index entry for the Chebyshev points now resolves to their definition, the Legendre points carry a defining footnote at first use, and the early exercise that reaches for the Chebyshev points before they are introduced now forward-references the definition.

    On one technical point the Referee's insight and the book's own argument converge more closely than the report allows. The Referee observes, correctly, that Computational Étude 4.1 in isolation studies only uniformly random nodes, and that the real difficulty is independent sampling, which lets points fall arbitrarily close together. That is precisely the conclusion the manuscript draws, in the very next étude. Computational Étude 4.2 samples random nodes carrying the correct arcsine marginal and finds that they perform even worse, concluding, in its own words, that "the success of Chebyshev nodes comes not from the arcsine marginal distribution, but from the deterministic structure that guarantees minimum node separation". I have adopted the Referee's further suggestion in full: the framing of Étude 4.1 now states explicitly that it is the first half of a two-part control that separates node density from node separation, and a new remark on determinantal point processes, which retain the arcsine marginal while adding exactly the repulsion the book identifies as essential, has been added, with a citation, to show that such processes confirm rather than contradict the book's thesis. I am grateful for a suggestion that let me sharpen the point in its modern language.

    The Referee's principal reservation concerns the book's engagement with the modern literature, and here I must respectfully distinguish what the book attempts from the standard against which it is measured. Each chapter closes with a section I have deliberately titled "A non-exhaustive literature overview", in all twenty-one instances; the adjective is chosen with care. These sections are not, and do not present themselves as, balanced or impact-weighted surveys of the state of the art. Their purpose is pedagogical: to show a graduate student, at the close of a chapter, that the technique just learned opens onto a living field, and to offer a few threads to pull. Read as the research surveys they explicitly decline to be, they will indeed seem uneven; read as what they are, they do the work intended. I have nonetheless taken the underlying concern seriously and checked afresh that the recent citations at issue are correctly and completely attributed, adding a few words of context where a reference had been left to stand too bare.

    Two specific claims I would gently correct, since they bear on the charge of imbalance. Fourier continuation, and the FFT-based treatment of non-periodic problems, is not absent: it is discussed in the chapter on Fourier grids and is the subject of a dedicated exercise that asks the reader to implement a least-squares or Gram-based Fourier extension and to compare it against the naive periodic derivative. And the ultraspherical sparse-spectral method of Olver and Townsend, far from being set on an equal footing with marginal work, is presented repeatedly and prominently as the paradigm shift it was, in the Afterword, in the linear-algebra appendix, and in the literature review of the chapter on coordinate transformations. I mention these not to dispute the Referee's good faith, but because the pattern of citation the report infers does not, I believe, survive a reading of the passages themselves.

    It may also help to set this report beside the others the series has received, since on the specific question of scholarship the readings genuinely diverge. The very features the fifth Referee finds wanting, the historical depth of the literature and the care behind it, were singled out for praise by two other referees who read the book in considerable detail: the first, whose expert report engaged the high-order and meshfree literature at length, and the second, whose six-page report commended the "historically deep literature surveys" by name. I record this not to set one referee against another, but because the six reports must be weighed together.

    On length, the Referee is quite right that the book is longer than any single course requires, and I have never supposed otherwise; roughly half of it makes a graduate course, and I know which half. The surplus is deliberate. The book is built to be modular, a resource from which an instructor assembles the course their own students need, whether that course leans toward fluid dynamics, toward eigenvalue problems, or toward time-dependent PDEs, rather than a single fixed syllabus; I have now made this intent explicit in the roadmap, which sets out the chapter dependencies, three course tracks, and the chapters that may be skipped. Much of the added length is the exercise set, which now numbers some three hundred and seventy problems, expanded, as the series editors will recall, at their own request as a condition of the book's acceptance into Texts in Applied Mathematics; the runnable code in three languages accounts for much of the rest, and is in any case archived in full in the companion repository.

    This brings me to the one characterisation I would respectfully resist. The Referee assesses the book as a reference "for graduate students or researchers" and finds it wanting in that role. It was never meant to fill it. The book is a graduate classroom textbook, and it says so in its own front matter: the Preface calls it "a book for the classroom and the laptop, not a monograph for the specialist's shelf", states in as many words that it is not "a research monograph, nor an advanced reference for specialists", and directs the specialist reader instead to Boyd, to Canuto and colleagues, to Shen, Tang and Wang, and to Trefethen for the complete theory and the state of the art; the Introduction repeats the point. That it has been conditionally accepted into a textbook series, Texts in Applied Mathematics, reflects the same understanding of what it is. Assessed against the standard of a research reference the book will fall short, as by design it must; assessed as the classroom text it declares itself to be, I believe it meets the standard that purpose sets.

    I am genuinely grateful to the Referee for the criticisms that have improved the book, and I have adopted every one that bears on its accuracy or its exposition. On the larger question of what kind of book this is, I ask only that it be read for what it sets out to do, and I am content to let the six reports, and the book itself, be weighed together.
  ],
)

#review-block(6,
  [
    #emph[The sixth Referee, writing after a complete reading, judges that "there is a very nice book lurking here" and recommends it once the text is shortened, its details polished, and the title changed. The report's substance is summarised here under the Referee's own headings.]

    #strong[Comparators and topic.] The book is read as a blend of Trefethen's mathematics, Boyd's discursive style, and the up-to-date Python/MATLAB/Julia approach of Driscoll and Braun. The Referee asks that the title be changed to make the subject, spectral methods for ODEs and PDEs, explicit, suggesting for instance "Spectral Methods: Computational Études".

    #strong[Length.] The central objection: at 857 pages, with some ninety pages of front matter (an eleven-page table of contents, a thirty-one-page list of figures, and further lists) before the Preface, the manuscript is far too long. The Referee would regard it as unpublishable beyond six hundred pages, would prefer five hundred, and points to occasional self-indulgence, citing a paragraph in the Afterword that recounts how that Afterword was drafted.

    #strong[Style and English.] Both are praised: the writing is "very good" and, as a discursive text, "more appealing than Boyd's" for being less eccentric; the English is "almost perfectly fluent". A few local lapses are noted: "obliged" after (4.15), "rigorousness" in §9.15, and the unexplained abbreviation "PNAS" in the same section.

    #strong[Wording.] Two openings overreach: Chapter 4 begins "Polynomial interpolation underlies every pseudospectral method", which is untrue for Fourier methods; and Chapter 6 reads as though every spectral approximation rested on a Fourier, rather than a Chebyshev, series.

    #strong[Figures.] The blue-toned palette is admired, but several early figures (Figures 2.2 to 2.4) are thought not compelling, and Figure 3.3 is faulted for drawing a "spectral basis" as a sawtooth when such functions are smooth.

    #strong[Code.] The three-language implementations, with all code available online, are called "excellent, reproducible science" that "students will love".

    #strong[Mathematics.] Mostly sound, with three imprecisions flagged: an Introduction passage describing the approximation of analytic functions as "explosive", faster than geometric, which does not hold in general; a suggestion in §4.2 that Runge's 1901 paper was about the function $1\/(1 + 25 x^2)$ rather than a general theory it merely illustrates; and a reference in §4.8.3 to a mean Lebesgue constant for random nodes, whose finiteness is questioned.

    #strong[History.] Judged good but not wholly reliable: the origin of spectral methods is told several times with a shifting cast (Lanczos, Kreiss and Oliger, Orszag, Gottlieb, Cooley and Tukey), the first telling omitting Orszag; the British pioneers (Clenshaw, Mason, Elliott, Fox) go unmentioned; the claim that the methods were prohibitive before the Cooley--Tukey FFT is doubted; and calling the ultraspherical method a "revolution" (§12.8) is thought excessive.

    #strong[Originality.] Considered limited, with several items close to Trefethen: Theorem 1 of §6.2 resembles a theorem of his spectral book, and Figures 12.1 and 12.3 resemble figures and an output of the same.

    #strong[Formatting and index.] The second paragraph of each chapter is not indented; and the index is thought superficial, most entries carrying a single page number, and is prefaced by a note that explains at too much length how it was built.

    In summary, the Referee sees "a very nice book lurking here" and asks the author to shorten the words, polish the details, and change the title.
  ],
  [
    I am warmly grateful to the sixth Referee for so encouraging a reading, and for the judgement that a very nice book is waiting to be brought out of the present manuscript. To be read in the company of Trefethen, Boyd, and Driscoll and Braun is a compliment I do not take lightly, and I have tried to earn it by acting on the report's advice throughout.

    On length, which is the heart of the report, I am in full agreement, and the manuscript has already been shortened substantially. The inline source-code listings, previously printed in all three languages for every étude, have been removed from the text and now live solely in the companion repository, where, as the Referee observes, the code is of most use; each étude points the reader to the exact files. The thirty-one-page List of Figures and the List of Tables have been removed, as has the Subject Index (on which more below), so that the reader reaches the Preface after a few pages rather than after ninety. The Afterword has been trimmed, and the self-referential passage recounting its own drafting, which the Referee justly singled out, has been deleted. Together these measures bring the manuscript from 857 pages down to 665, a reduction of some one hundred and ninety pages. I have deliberately preserved the mathematics, the computational experiments, the exercises, and the end-of-chapter literature surveys, which two of the other referees praised warmly; should the series editors wish to bring the book below six hundred pages, those optional surveys are the natural place to trim further, and I would be glad to do so.

    On the title, I am happy to make the spectral-methods subject explicit, and I leave the final wording gladly to you and to the series editors.

    The two chapter openings the Referee flagged have been corrected. Chapter 4 no longer claims that polynomial interpolation underlies every pseudospectral method, but only the methods for non-periodic problems, noting that their periodic counterparts rest on trigonometric interpolation; and Chapter 6 now develops the argument from smoothness to coefficient decay for Chebyshev as well as Fourier coefficients, since the substitution $x = cos theta$ carries every estimate from one to the other. The smaller lapses are also mended: "rigorousness" is now "rigour", "PNAS" is spelled out at first use as the Proceedings of the National Academy of Sciences, and the sentence after (4.15) that used "obliged", and read as a contradiction, has been rewritten.

    Figure 3.3 has been redrawn so that the spectral basis function appears as a smooth curve, a Chebyshev polynomial, rather than a sawtooth. I am grateful, too, for the kind words on the figures' blue palette, and I will give the earliest plots the further attention the Referee's remark invites.

    The three mathematical imprecisions are corrected. The Introduction now separates the geometric convergence of functions analytic on a neighbourhood of the interval from the faster-than-geometric convergence reserved for entire functions, and no longer calls the analytic case "explosive". Section 4.2 now presents Runge's 1901 work as a general study of equispaced interpolation, illustrated by the famous example $1\/(1 + 25 x^2)$ rather than being about it. And §4.8.3 now leads with the median Lebesgue constant and states plainly that, because independent sampling can place two nodes arbitrarily close, the population mean may fail to be finite; the sample mean is offered only as a diagnostic of that heavy tail. The Referee's instinct on this last point was exactly right.

    On the history, the origin story has been made consistent, and Orszag, whose part in making the methods known was indeed decisive, now appears in the first telling as well as the later ones. The British pioneers, Clenshaw and, after him, Fox, Mason, and Elliott, are credited for turning Chebyshev series into a practical instrument, with a reference to Fox and Parker. The suggestion that the methods were prohibitive before the Cooley--Tukey FFT has been softened: for the modest problem sizes of the 1960s the direct transform was manageable, and it is at large scale that the FFT proved decisive. And the ultraspherical method is no longer described as a "revolution".

    On originality, the debt to Trefethen is real and, I hope, openly owned: the book is in part an homage, and it re-presents the classical results and canonical examples of the field for a new generation and a new computational setting. Where that debt was insufficiently marked I have marked it. Theorem 1 of §6.2, a classical Paley--Wiener result, now carries an explicit citation to the formulation in Trefethen's spectral book; and the captions of Figures 12.1 and 12.3 now state that they follow, respectively, his polar-grid illustration and his Program 28. What the book offers of its own is not new theorems, but a unified and fully reproducible treatment in three languages, the étude format, and a large graded set of exercises.

    Finally, the paragraph-indentation fault has been corrected in the template, and the index, which the Referee rightly found superficial, has been removed for this revision, together with the over-informative note that headed it; a professional index can be prepared at the page-proof stage, as is customary for the series.

    I am grateful to the Referee for a report that has made the book shorter, sharper, and more accurate, and for the encouragement to bring out the better book that, as the report kindly puts it, lies within this one.
  ],
)
