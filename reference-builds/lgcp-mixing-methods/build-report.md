# Bygg- och verifieringsrapport

Datum: 2026-11-10. Status: färdig riktad litteraturöversikt och metodrekommendation; ingen ny samplarimplementation eller empirisk algoritmjämförelse.

## Leverans

* [Svensk rapport](lgcp_methods_review.pdf), sju sidor med prioritering, modellkoppling, aGrad-z-formler, CV-diskussion och 14 referenser.
* [TeX-källa](lgcp_methods_review.tex).
* [Källguide med inbäddad BibTeX](reference.md), inklusive samtliga filvägar, sidlokaliseringar, versioner och kontrollsummor.
* [Maskinläsbar filförteckning](source-inventory.json).
* 15 bevarade PDF-filer för 14 arbeten i `papers/`, med 15 råextraktioner och en separat OCR i `source-text/`.

Källnycklar: `murray2010elliptical`, `knorrheld2002block`, `titsias2018auxiliary`, `cotter2013functions`, `beskos2011hilbert`, `girolami2011riemann`, `cui2016dili`, `papaspiliopoulos2007parametrization`, `syed2022tempering`, `rue2009inla`, `andrieu2009pseudo`, `vehtari2021rhat`, `vehtari2024psis`, `vehtari2017loo`.

## Sökning och urval

Riktad sökning efter metoder för beroende Gaussiska fält, LGCP, gles precision, blockuppdatering, gradientbaserad MCMC, parametrisering, temperering och MCMC-baserad CV. Förlags-, konferens-, författar- och institutionskällor samt arXiv-fulltexter användes. Publicerad metadata verifierades mot primära bibliografiska poster; metodpåståenden kontrollerades i fulltext. Detta är inte en systematisk review och kan inte användas som belägg för att ingen senare eller bättre metod finns.

Särskilt viktiga primära metadataposter är PMLR:s sida för Murray m.fl., Wiley/OUP:s tidskriftsposter för Knorr-Held–Rue och Titsias–Papaspiliopoulos, Elseviers artikelpost för Cui m.fl., UCL:s post för Beskos m.fl., JMLR:s post för PSIS och Aaltos institutionspost för LOO-artikeln. DOI/stabila länkar finns i BibTeX och hämtadresser i filförteckningen.

## Extraktion och versioner

Poppler `pdftotext` 26.04.0 med `-layout -enc UTF-8`; ettbaserade PDF-sidpositioner och form feed behålls. Alla sparade PDF- och textkontrollsummor verifierades efter arkivering. Alla intakta källor har lika många textsidsegment som PDF-sidor. Detta visar täckning av sidgränser, inte att varje tecken är korrekt.

Två läsbarhetsproblem hanterades:

1. Knorr-Held–Rue: råtexten hade skadad teckenkodning. Samtliga 39 sidor OCR-behandlades separat med Tesseract 5.5.2 efter rendering i skala 2.5. Matematiska uttryck ska fortfarande hämtas ur PDF-bilden.
2. Girolami–Calderhead: första kopian var skadad, med 36 textsegment trots 92 deklarerade sidor. En fullständig 92-sidig publicerad kopia hämtades och används för belägg. Den skadade kopian är tydligt märkt och bevarad för proveniens.

Syed-filen har kraftigt upprepad figurtext men full sidsegmenttäckning; relevanta prosasidor granskades visuellt. Skillnader mellan förtrycksår och publiceringsår, reprint- och tidskriftssidor samt diskussionsbilagor är redovisade i referensguidens avsnitt 0. Särskilt verifierades DOI `10.1111/rssb.12269`, Cui-artikelns sidor 109–137 och PSIS-artikelns sidor 1–58. Ingen likställning av olika versionsnumreringar har gjorts.

## Innehållskontroll

Visuellt kontrollerade källsidor (ettbaserade PDF-positioner):

| Källa | Sidor och användning |
|---|---|
| Murray m.fl. | 3: figur 2 och beräkningskostnad |
| Knorr-Held–Rue | 3, 6: resultatens empiriska omfattning och blockkonstruktion |
| Titsias–Papaspiliopoulos | 10, 17–18, 39: definition av aGrad-z, LGCP-upplägg/tabell och algoritm 1 |
| Cotter m.fl. | 11–12: pCN, acceptans och pCNL |
| Beskos m.fl. | 14: HMC-konstruktionen i tabell 2 |
| Girolami–Calderhead | 23: LGCP-täthet och konstant fältmetrik |
| Cui m.fl. | 11, 19: informativt delrum och förslag |
| Papaspiliopoulos m.fl. | 5: icke-centrering, definition och informationsberoende |
| Syed m.fl. | 2: temperaturväg och utvidgat mål |
| Rue m.fl. | 3, 23: latent Gaussisk struktur och LGCP-faktorisering |
| Andrieu–Roberts | 6: utvidgat mål och marginalegenskap |
| Vehtari m.fl. 2021 | 4: rekommendationer om R-hat, kedjor och ESS |
| Vehtari m.fl. 2024 | 5: PSIS-algoritm och relativ MCMC-effektivitet |
| Vehtari m.fl. 2017 | 8: datadelning och gruppvis CV |

Kompletterande sammanhang lästes i den bevarade texten enligt lokaliseringarna i referensguiden. Interna definitioner för de återgivna algoritmerna följdes. Externa resultat som källorna hänvisar till men som inte granskats separat identifieras där; deras bevis eller egna experiment förs inte vidare som självständigt verifierade resultat.

Kodbeskrivningen kontrollerades mot `R/lgcp.R`, `R/lgcp_sampler.R` och `R/lgcp_blocks.R`. Blandningsresultaten kommer från `inst/reports/lgcp/field_mixing_interpretation.tex` och tillhörande diagnostik, inte från nya körningar. Rapporten skiljer uttryckligen mellan observerade och latenta antal, mellan korrekt MH och approximation samt mellan betingad och gemensam blandning.

Definitioner och notation kontrollerades i en egen genomläsning: x, P, B, a, f och mu definieras före gradientalgoritmen; f innehåller inte priorn. Aktivitetsmarginaliseringen och den prediktiva estimanden märks som direkt algebra för vår modell. Ingen oberoende matematisk bevisgranskning har gjorts, och ingen generisk funktionsrumssats intygas för river-modellen.

## Teknisk och visuell kontroll

* 14 unika BibTeX-poster extraherades från Markdown till byggkatalogen och accepterades av BibTeX 0.99d med `plainnat`.
* Samtliga 14 citerade nycklar matchar bibliografin. Inga absoluta lokala filfält förekommer i BibTeX.
* Rapporten byggdes med pdfLaTeX; sista bygget saknar odefinierade referenser, varningar och överfulla/underfulla rader.
* Samtliga sju rapportsidor renderades och granskades visuellt: tabell, matematik, sidnummer, hänvisningar och referenslista är läsbara utan klippning eller överlappning.
* Lokala dokumentlänkar och PDF-/textkontrollsummor kontrollerades.
* Det fristående byggskriptet kördes. Slutfilens text på alla sju sidor överensstämmer med den visuellt granskade byggversionen. PDF SHA-256: `fd6f436fcc55f3b48733681d745c5a9c3b7302d7903d7b24d9e973ad06ad0c96`.
* `.Rbuildignore` fick en avgränsad rad för `reference-builds`, så att källarkivet inte följer med i R-paketets byggfil. Inferenskoden ändrades inte av litteraturarbetet.

## Återskapa PDF

Kör `python reference-builds/lgcp-mixing-methods/build.py` från projektroten. Det kräver `pdflatex` och `bibtex`. Bibliografin läses från det inbäddade blocket i `reference.md`; temporära byggfiler rensas automatiskt. Ingen separat `.bib`-leverans behövs.

## Begränsningar och nästa användning

Redo som underlag för metodval och experimentdesign. Inte redo som stöd för en påstådd uppmätt hastighetsvinst eller konvergensgaranti. Samtliga nya algoritmers prestanda i den fullständiga modellen återstår att testa. Inte heller har INLA, variationsinferens eller pseudo-marginal integration implementerats för den utökade modellen. Den rekommenderade ordningen och försöksplanen är författarens tillämpningsbedömning, separat från de citerade arbetenas egna resultat.
