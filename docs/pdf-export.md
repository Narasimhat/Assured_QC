# PDF export of the figures

Every figure of Assured QC can be saved as a vector PDF as well as SVG and PNG. The PDF is made by `src/svgPdf.js` from the same SVG that is shown on screen, so the two always agree; it needs no extra package and no server (nothing leaves the browser).

## Where

- **Browser, run results:** the buttons under each figure: Clone selection, Base-call alignment, Chromatograms (all clones aligned), Summary graph, Protein products, each with an (SVG), (PNG) and (PDF) button. For the selected sample: Aligned chromatogram and Protein products in the detail view, and Figure (PDF) above it.
- **Command line:** `run` writes `summary.pdf`, `clone_selection.pdf`, `clones_alignment.pdf`, `clones_aligned.pdf` (and `protein_products.pdf` when there is a gene record) next to the SVGs, and one PDF per clone in `chromatograms/` with `--clone-figures`. `--no-pdf` (put it last) turns this off. `figure --out figure.pdf` writes the single-sample figure as a PDF.

## What the PDF contains

- One page, as large as the figure (one SVG pixel is 0.75 pt, so a 980 px wide figure is 735 pt, 10.2 inches, wide). Everything is vector, so it scales without loss and can be edited in Illustrator or Inkscape; chromatogram traces are polylines, not pictures.
- Text is set in the standard Helvetica and Courier faces, which are not embedded. Every PDF viewer has them and they print and edit normally; if a journal asks for embedded fonts, open the PDF in Illustrator or Acrobat and save with fonts embedded, or use the SVG.
- The minus sign, en dash, ellipsis, middle dot and non-breaking space are drawn as such. A character that has no glyph in these fonts is replaced by a look-alike (for example `>=` for the greater-or-equal sign) or by `?`; sample names in other scripts would be affected.

## What is supported

The converter handles the elements the figures of the app use: `svg`, `g` with `transform`, `rect` (also rounded), `circle`, `ellipse`, `line`, `polyline`, `polygon`, `path` (all commands, including arcs), `text` (size, bold, italic, anchor), fill and stroke colour, opacity, stroke width, dashes, line joins and caps, and the tiling pattern of the protein figure. Any other element (an embedded image, a filter, a gradient) is skipped and reported by `convertSvg` in `unsupported`; the tests fail if a figure of the app contains one.

## How it was checked

Each figure of the app is converted from simulated traces in `test/svg-pdf.test.js` (page size, one text operator per drawn text, cross-reference table, no unsupported element, no lost character, compressed and plain page content identical). In addition the PDFs were rendered and laid over the SVG renders of the same figures (clone selection, summary, aligned chromatograms, base-call alignment, protein products with the hatch fill); the layouts agree and the remaining differences are anti-aliasing.
