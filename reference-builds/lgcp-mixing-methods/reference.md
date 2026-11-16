# Metoder för bättre blandning i river-LGCP

Riktad litteraturöversikt, 2026-11-10. Den sammanhängande svenska rapporten finns i [PDF](lgcp_methods_review.pdf), med redigerbar [TeX-källa](lgcp_methods_review.tex). Här finns källkartan och underlaget för återanvändning. Detta är ett urval av 14 primärkällor, inte en systematisk översikt över all LGCP-litteratur eller alla senare MCMC-metoder.

**Bedömning för odispers:** jämför först befintliga globala ESS- och Laplace-MH-steg. Prioritera sedan aGrad-z med gles precision och förkonditionerad HMC. Utvärdera både den betingade fältfördelningen vid fixerade antal och hela samplaren med latenta antal. INLA och PSIS löser andra delproblem och ersätter inte denna blandningskontroll.

## Källkarta

| Fråga | Primärkälla och kontrollerad lokalisering | Avsnitt nedan |
|---|---|---|
| Varför större block? | Knorr-Held–Rue, reviderat förtryck: abstract, PDF s. 3; §2, s. 6–7 | 2 |
| Vad gör elliptisk slice sampling? | Murray m.fl.: figur 2 och §2.2, PDF s. 3 | 2 |
| Vilken ny gradientmetod är mest direkt relevant? | Titsias–Papaspiliopoulos: §3.3, s. 8–11; §5.2/tabell 2, s. 17–18; supplement algoritm 1, s. 39 | 3 |
| Prioranpassning och HMC | Cotter m.fl.: §4.2–4.3, s. 10–12; Beskos m.fl.: §3.3/tabell 2, s. 12–14; Girolami–Calderhead: §9, s. 23–24 | 4 |
| Anpassning till informativa riktningar | Cui m.fl.: §3.2, s. 11; §3.3.1, s. 19 | 5 |
| Koppling mellan latenta variabler | Papaspiliopoulos m.fl.: §2 (figur 3–4), s. 5 | 6 |
| Numerisk marginalisering | Andrieu–Roberts: §2, ekvation (2.2)–(2.4), s. 6 | 6 |
| Barriärer mellan regimer | Syed m.fl.: §1.1, s. 2 | 7 |
| Kan INLA ersätta MCMC? | Rue m.fl.: §1.3, s. 3; LGCP-exempel §5.5 börjar s. 23 | 8 |
| Blandning och pålitlig CV | Vehtari m.fl. 2021: §2, s. 4; 2017: §2.4, s. 8; 2024: algoritm 1/ekv. (6), s. 5 | 9 |

Alla sidtal ovan är **ettbaserade PDF-positioner i de sparade versionerna**, inte tidskriftens tryckta sidnummer. PDF-filerna är auktoritativa för formler. Sökbar text är ett hjälpmedel. Källkontroll av metodbeskrivningar ska inte förväxlas med oberoende bevisgranskning.

## 0. Källfiler, versioner och bibliografi

Publicerad metadata används i BibTeX. Exakta sidlokaliseringar avser de fulltexter som anges nedan; likalydande avsnittsnummer i en annan version ska kontrolleras på nytt.

* Knorr-Held–Rue: publicerad 2002 i *Scandinavian Journal of Statistics*, 29(4), 597–614. Sparad fulltext är rapport 210, första version juli 2000, reviderad maj 2001, inklusive institutionsomslag. Tidskrifts-DOI är `10.1111/1467-9469.00308`; repository-DOI `10.5282/ubm/epub.1600` identifierar rapportkopian. Den råa textextraktionen är gravt teckenkodningsskadad. En separat OCR av samtliga 39 sidor är sökbar men opålitlig för matematiska symboler.
* Titsias–Papaspiliopoulos: publicerad 2018; läst arXiv:1610.09641v3, 25 januari 2018, 41 sidor inklusive supplement. DOI slutar **12269**. Supplementets algoritm 1 är den kontrollerade källan för aGrad-z-formlerna i rapporten.
* Cotter m.fl.: sparad IMS-reprint med ett extra Manchesteromslag. PDF-sida 11 är reprintens sida 10. Bibliografin använder tidskriftens sidor 424–446.
* Beskos m.fl.: författarversion daterad 5 maj 2011; bibliografin använder publicerade uppgifter 2011, 121(10), 2201–2230.
* Girolami–Calderhead: den första hämtade Cornellkopian deklarerar 92 sidor men gav endast 36 textsidsegment och flera läsfel. Den bevaras som en skadad källa och används **inte** för sidlokalisering eller innehållsbelägg. Den senare författar-/diskutantkopian från Julien Cornebise innehåller samtliga 92 sidor och används här. PDF-sida 23 motsvarar tryckt sida 145. Använd inte arXiv:1001.4460 som identifierare för denna artikel; den tidigare arXiv:0907.1100-versionen har annan titel och författarlista och ingår inte i paketet.
* Cui–Law–Marzouk: läst arXiv:1411.3688v2, 17 oktober 2015; publicerad 2016. Sidintervallet är **109–137**. Författarordningen är Cui, Law, Marzouk.
* Papaspiliopoulos–Roberts–Sköld och Andrieu–Roberts: IMS-reprints med egen sidnumrering. Publicerade sidintervall anges i bibliografin.
* Syed m.fl.: läst arXiv:1905.02939v4, juli 2021, 74 sidor inklusive appendicer. Tidskriftsnumret är från **2022**, även om tidig nätpublicering/repository anger 2021. Figurtext är mycket repetitiv i råextraktionen; matematiska/prosaiska passager som används har kontrollerats visuellt.
* Rue–Martino–Chopin: läst 40-sidigt förtryck från 2008; publicerad artikel med diskussion är från 2009, s. 319–392. LGCP-avsnittets nummer och sidpositioner avser förtrycket.
* Vehtari m.fl. 2021: läst arXiv:1903.08008v5, 22 juni 2021, 26 sidor; tidskriftens metadata inkluderar diskussion och sidorna 667–718.
* Vehtari–Gelman–Gabry: läst förtryck daterat 1 september 2016, 28 sidor; publicerad årgång 2017. DOI innehåller 2016.
* Vehtari m.fl. 2024: JMLR:s publicerade PDF, 25(72), **1–58**. En äldre uppgift om 57 sidor ska inte användas för denna version.

### BibTeX

```bibtex
@inproceedings{murray2010elliptical,
  author = {Murray, Iain and Adams, Ryan and MacKay, David},
  title = {Elliptical slice sampling},
  booktitle = {Proceedings of the Thirteenth International Conference on Artificial Intelligence and Statistics},
  series = {Proceedings of Machine Learning Research},
  volume = {9}, pages = {541--548}, year = {2010}, publisher = {PMLR},
  url = {https://proceedings.mlr.press/v9/murray10a.html}
}
@article{knorrheld2002block,
  author = {Knorr-Held, Leonhard and Rue, H{\aa}vard},
  title = {On Block Updating in {Markov} Random Field Models for Disease Mapping},
  journal = {Scandinavian Journal of Statistics},
  volume = {29}, number = {4}, pages = {597--614}, year = {2002},
  doi = {10.1111/1467-9469.00308}, url = {https://doi.org/10.1111/1467-9469.00308}
}
@article{titsias2018auxiliary,
  author = {Titsias, Michalis K. and Papaspiliopoulos, Omiros},
  title = {Auxiliary gradient-based sampling algorithms},
  journal = {Journal of the Royal Statistical Society: Series B},
  volume = {80}, number = {4}, pages = {749--767}, year = {2018},
  doi = {10.1111/rssb.12269}, url = {https://doi.org/10.1111/rssb.12269}
}
@article{cotter2013functions,
  author = {Cotter, S. L. and Roberts, G. O. and Stuart, A. M. and White, D.},
  title = {{MCMC} Methods for Functions: Modifying Old Algorithms to Make Them Faster},
  journal = {Statistical Science},
  volume = {28}, number = {3}, pages = {424--446}, year = {2013},
  doi = {10.1214/13-STS421}, url = {https://doi.org/10.1214/13-STS421}
}
@article{beskos2011hilbert,
  author = {Beskos, A. and Pinski, F. J. and Sanz-Serna, J. M. and Stuart, A. M.},
  title = {Hybrid {Monte Carlo} on {Hilbert} spaces},
  journal = {Stochastic Processes and their Applications},
  volume = {121}, number = {10}, pages = {2201--2230}, year = {2011},
  doi = {10.1016/j.spa.2011.06.003}, url = {https://doi.org/10.1016/j.spa.2011.06.003}
}
@article{girolami2011riemann,
  author = {Girolami, Mark and Calderhead, Ben},
  title = {Riemann manifold {Langevin} and {Hamiltonian Monte Carlo} methods},
  journal = {Journal of the Royal Statistical Society: Series B},
  volume = {73}, number = {2}, pages = {123--214}, year = {2011},
  doi = {10.1111/j.1467-9868.2010.00765.x}, url = {https://doi.org/10.1111/j.1467-9868.2010.00765.x}
}
@article{cui2016dili,
  author = {Cui, Tiangang and Law, Kody J. H. and Marzouk, Youssef M.},
  title = {Dimension-independent likelihood-informed {MCMC}},
  journal = {Journal of Computational Physics},
  volume = {304}, pages = {109--137}, year = {2016},
  doi = {10.1016/j.jcp.2015.10.008}, url = {https://doi.org/10.1016/j.jcp.2015.10.008}
}
@article{papaspiliopoulos2007parametrization,
  author = {Papaspiliopoulos, Omiros and Roberts, Gareth O. and Sk{\"o}ld, Martin},
  title = {A General Framework for the Parametrization of Hierarchical Models},
  journal = {Statistical Science},
  volume = {22}, number = {1}, pages = {59--73}, year = {2007},
  doi = {10.1214/088342307000000014}, url = {https://doi.org/10.1214/088342307000000014}
}
@article{syed2022tempering,
  author = {Syed, Saifuddin and Bouchard-C{\^o}t{\'e}, Alexandre and Deligiannidis, George and Doucet, Arnaud},
  title = {Non-Reversible Parallel Tempering: A Scalable Highly Parallel {MCMC} Scheme},
  journal = {Journal of the Royal Statistical Society: Series B},
  volume = {84}, number = {2}, pages = {321--350}, year = {2022},
  doi = {10.1111/rssb.12464}, url = {https://doi.org/10.1111/rssb.12464}
}
@article{rue2009inla,
  author = {Rue, H{\aa}vard and Martino, Sara and Chopin, Nicolas},
  title = {Approximate {Bayesian} inference for latent {Gaussian} models by using integrated nested {Laplace} approximations},
  journal = {Journal of the Royal Statistical Society: Series B},
  volume = {71}, number = {2}, pages = {319--392}, year = {2009},
  doi = {10.1111/j.1467-9868.2008.00700.x}, url = {https://doi.org/10.1111/j.1467-9868.2008.00700.x}
}
@article{andrieu2009pseudo,
  author = {Andrieu, Christophe and Roberts, Gareth O.},
  title = {The pseudo-marginal approach for efficient {Monte Carlo} computations},
  journal = {The Annals of Statistics},
  volume = {37}, number = {2}, pages = {697--725}, year = {2009},
  doi = {10.1214/07-AOS574}, url = {https://doi.org/10.1214/07-AOS574}
}
@article{vehtari2021rhat,
  author = {Vehtari, Aki and Gelman, Andrew and Simpson, Daniel and Carpenter, Bob and B{\"u}rkner, Paul-Christian},
  title = {Rank-normalization, folding, and localization: An improved {$\widehat{R}$} for assessing convergence of {MCMC} (with discussion)},
  journal = {Bayesian Analysis},
  volume = {16}, number = {2}, pages = {667--718}, year = {2021},
  doi = {10.1214/20-BA1221}, url = {https://doi.org/10.1214/20-BA1221}
}
@article{vehtari2024psis,
  author = {Vehtari, Aki and Simpson, Daniel and Gelman, Andrew and Yao, Yuling and Gabry, Jonah},
  title = {Pareto Smoothed Importance Sampling},
  journal = {Journal of Machine Learning Research},
  volume = {25}, number = {72}, pages = {1--58}, year = {2024},
  url = {https://www.jmlr.org/papers/v25/19-556.html}
}
@article{vehtari2017loo,
  author = {Vehtari, Aki and Gelman, Andrew and Gabry, Jonah},
  title = {Practical {Bayesian} model evaluation using leave-one-out cross-validation and {WAIC}},
  journal = {Statistics and Computing},
  volume = {27}, number = {5}, pages = {1413--1432}, year = {2017},
  doi = {10.1007/s11222-016-9696-4}, url = {https://doi.org/10.1007/s11222-016-9696-4}
}
```

## 1. Modell och tillämpningsgräns

Använd denna referens för val av fältuppdatering, frågor om lång autokorrelation, jämförelse av beräkningstid mot effektiv stickprovsstorlek och bedömning av CV under MCMC-beroende.

Kodkontroll: `R/lgcp_sampler.R`, funktionerna `lgcp_poisson_loglik`, `lgcp_field_log_density`, `lgcp_field_laplace` och `lgcp_field_step`; dispatch i `R/lgcp.R`; lokala block i `R/lgcp_blocks.R`. Den gemensamma kontinuerliga variabeln består av intensitetskoefficienter och fält. Villkorat på punktantalen är dess täthet proportionell mot exp(f(x)) N(x;0,P^-1), där f är Poissonloglikelihooden. Gradient och negativ Hessian i PDF-rapporten följer direkt från denna täthet. Antal och märken är diskreta och kräver separat behandling.

Den huvudsakliga tillämpningsgränsen är **observerade kontra latenta punktantal**. LGCP-fältresultat kan vara direkt relevanta för vår betingade uppdatering utan att lösa beroendet i hela den hierarkiska modellen. Fältets precision är fixerad inom en körning. Skalan/räckvidden behöver därför inte samplas tillsammans med fältet i den nuvarande jämförelsen.

## 2. Stora block och elliptisk slice sampling

**Källstatus:** metodkonstruktion och empiriska jämförelser. **Kontroll:** Knorr-Held–Rue abstract PDF s. 3 och §2 s. 6–7; Murray m.fl. figur 2/§2.2 PDF s. 3. Ingen fristående granskning av invariansbevis eller av Rue (2001), som blockartikeln hänvisar till för Gaussisk simulering. Detaljer ur det externa beviset återges därför inte.

Knorr-Held–Rue rapporterar stora vinster av gemensamma block i sina sjukdomskarteringsexempel, särskilt när fält och hyperparameter uppdateras tillsammans. Murray m.fl. ger en algoritm som söker längs en Gaussisk priorellips med en likelihoodbaserad slicetröskel. Att undvika slutlig förkastning betyder inte att successiva drag blir oberoende.

**Vår bedömning:** benchmarka de redan befintliga globala alternativen före större nyutveckling. Små lokala block kan bidra till långsam rörelse i rumsliga medelvärden. Det är ännu ingen isolerad orsaksanalys. Laplace-MH i koden har acceptanskorrigering; att ersätta den med okorrigerade Gaussiska drag skulle ändra inferensen.

## 3. Auxiliary gradient sampling

**Källstatus:** MCMC-konstruktion och experiment. **Kontroll:** Titsias–Papaspiliopoulos §3.3 PDF s. 8–11, §5.2 s. 17–18, supplement algoritm 1 s. 39. Figur-/algoritmsidorna har granskats visuellt. Algoritmens notation och Gaussiska priorantagande följs tillbaka till §3.3; ingen generell konvergenssats återges.

Deras LGCP-experiment använder 64×64 celler, observerade antal, fixerade hyperparametrar och tio upprepningar. Minsta ESS/s i tabell 2 är 0.40 för aGrad-z, 1.35 för mGrad och 0.08 för Ellipt. Detta stöder ett metodval för test, inte en förutsägelse om en bestämd hastighetsvinst här.

**Vår bedömning:** aGrad-z kan uttryckas genom faktorisering av P+2I/δ; algoritmen i PDF-rapporten är supplementets algoritm 1 med C^-1=P. Inga täta kovariansmatriser behöver bildas för just den implementationen. Fast P och δ möjliggör återanvändning av faktorisering, medan gradienten ändras när antal och fält ändras. Skalning, acceptans och faktisk faktorfyllnad måste mätas. Denna glesa anpassning har inte implementerats eller benchmarkats i översikten.

## 4. Prioranpassade förslag och HMC

**Källstatus:** algoritmer med funktionsrumsanalys respektive LGCP-experiment. **Kontroll:** Cotter m.fl. §3 PDF s. 8 och §4.2–4.3 s. 10–12, särskilt (4.8) s. 11; Beskos m.fl. §3.3/tabell 2 s. 12–14; Girolami–Calderhead §9 s. 23–24 i den fullständiga publicerade kopian.

pCN bevarar den Gaussiska priorn i förslaget och accepterar med likelihoodkvoten. pCNL använder även likelihoodgradienten. Hilbert-HMC bygger vidare på en Gaussisk referensfördelning och Hamiltonsk dynamik. Girolami–Calderheads LGCP-fält har en konstant metrik vid fixerade hyperparametrar. Deras §9 återanvänder ett exempel från Christensen m.fl. (2005); den ursprungliga artikeln har inte granskats separat och dess egna resultat åberopas inte här.

**Vår bedömning:** priorvitning eller en fast krökningsanpassning är rimliga HMC-piloter för (theta,g)|N. Pappersresultat om välformulerade funktionsrumsmetoder är inte generella löften om låg autokorrelation; inga regularitetsvillkor eller dimensionsoberoende satser förs över till vår river-modell utan ny kontroll. Vanlig HMC kan inte direkt uppdatera diskreta antal.

## 5. Likelihood-informed subspace / DILI

**Källstatus:** konstruktion och experiment, inte en påvisad egenskap hos vår modell. **Kontroll:** Cui–Law–Marzouk §3.2 PDF s. 11 och §3.3.1 s. 19. Definitionen av det informativa delrummet följs till den priorförkonditionerade krökningen; inga optimala lågrangssatser från artikelns externa referenser återges.

DILI använder olika förslag i informativa riktningar och deras komplement. **Vår bedömning:** undersök om ett litet delrum faktiskt fångar krökning och långsamma funktioner. Antalet observerade byar begränsar inte utan vidare rangen hos den betingade Poissoninformationen när latenta cellantal behandlas som givna. En billig lågdimensionell korrigering är en hypotes att testa.

## 6. Omparametrisering, kollapsning och pseudo-marginal MCMC

**Källstatus:** generella konstruktioner och illustrerande modelljämförelser. **Kontroll:** Papaspiliopoulos–Roberts–Sköld §2 (figur 3–4) PDF s. 5; Andrieu–Roberts §2, ekvation (2.2)–(2.4), PDF s. 6. I det senare fallet kontrolleras den utvidgade täthetens marginalegenskap och kravet på en korrekt vikt-/likelihoodskattare; inga ergodicitetssatser används.

Icke-centrering kan minska hierarkiskt beroende, men är inte universellt bättre. Med fast Q är vår vitning en linjär koordinattransformation. Pseudo-marginalmetoden använder en icke-negativ, väntevärdesriktig skattare för att bevara rätt marginalmål i ett utvidgat rum; det räcker inte med en godtycklig noisy likelihood.

**Vår bedömning:** ett gemensamt fält–antalsförslag kräver korrekt MH-kvot för hela berörda mål- och förslagsfördelningen. Aktivitetsindikatorn kan däremot summeras bort exakt med den tvåfallsformel som ges i rapporten. Den formeln är vår direkta algebra, inte ett resultat hämtat ur dessa artiklar. Stabil numerisk integration över latenta platsmönster återstår att utveckla. Förväntad exponering är inte likvärdig med att integrera den ursprungliga bylikelihooden.

## 7. Parallel tempering

**Källstatus:** metodkonstruktion. **Kontroll:** Syed m.fl. §1.1 PDF s. 2, där temperaturväg, produktmål och lokala kedjor definieras. Inga asymptotiska effektivitetsgränser eller garantier för vår modell återges.

Metoden utbyter tillstånd mellan tempererade fördelningar för att underlätta passage mellan regimer. **Vår bedömning:** relevant först om barriärer kvarstår efter bättre kontinuerliga steg. Olika bredder i okonvergerade kedjor bevisar inte multimodalitet. Valet av temperaturväg, stöd vid likelihood noll, antal temperaturer och kostnad per fold måste bestämmas och valideras.

## 8. INLA och approximation

**Källstatus:** approximativ inferens för latenta Gaussiska modeller. **Kontroll:** Rue–Martino–Chopin §1.3 PDF s. 3; §5.5 från s. 23 visar artikelns LGCP-tillämpning. De antaganden som används här är Gaussiskt latent fält och observationsfaktorisering; inga numeriska noggrannhetsresultat överförs.

**Vår bedömning:** vår hela modell är inte direkt den vanliga observerade-LGCP-modellen eftersom platsantal och märken är latenta och flera byar delar exponering från samma platser. En betingad Laplaceapproximation kan vara användbar. Huruvida en större INLA- eller variationsbaserad omformulering blir lämplig är öppet; variationslitteraturen har inte särskilt granskats i detta urval.

## 9. Blandning, prediktiva integraler och CV

**Källstatus:** diagnostiska rekommendationer och importance-samplingmetoder. **Kontroll:** Vehtari m.fl. 2021 §2 PDF s. 4; Vehtari–Gelman–Gabry 2017 §2.4 s. 8; Vehtari m.fl. 2024 algoritm 1 samt ekv. (6) och efterföljande definition av relativ MCMC-effektivitet s. 5.

Ranknormaliserad R-hat under 1.01, flera kedjor och funktionsanpassad ESS/MCSE är användbara kontroller, inte bevis för full utforskning. Gruppvis utelämning ställer större krav på importancevikterna än punktvis LOO. PSIS hanterar viktsvansar men kan inte återställa posteriora regioner som samplaren aldrig har besökt; dess MCSE behöver ta hänsyn till MCMC-beroende.

**Vår bedömning:** följ fältets stora skalor och de verkliga prediktiva integralerna i varje fold. Logaritmen av den posteriora medellikelihooden är inte samma estimand som medelvärdet av loglikelihooden. Använd stabila kedjor och Monte Carlo-fel för den faktiska aggregerade poängen; skilj detta från variation mellan datamängder/foldar. Behåll separat foldanpassning som referens för den geografiska CV-uppgiften.

## 10. Kontroller, öppna frågor och återanvändning

Översikten stöder metodval och design av ett jämförande experiment. Den visar inte att någon alternativ samplare redan är snabbare eller konvergerad i odispers. Huvudfrågorna är om långsamheten domineras av fältets rumsliga skalor, av fält–antalskoppling eller av barriärer i den gemensamma posteriorn, och hur det påverkar de prediktiva funktionerna.

Kontrollerad täckning: listade metodkonstruktioner, förutsättningar för den direkta betingade tillämpningen, utvalda originalexperiment samt diagnostisk tolkning. Ej granskad täckning: fullständiga bevis, fullständig diskussionslitteratur, alla senare samplare, INLA-implementation för den utökade modellen eller variationsinferens. En oberoende matematisk bevisrevision har inte gjorts.

Fullständig filförteckning och kontrollsummor följer nedan. Läs även [bygg- och verifieringsrapporten](build-report.md). För sökning, från denna katalog:

```sh
rg -n -i 'aGrad-z|Log-Gaussian|likelihood-informed' source-text/*.txt
```

För Knorr-Held–Rue bör OCR-filen användas för ordsökning och original-PDF:en för exakta uttryck. Använd den fullständiga `girolami-calderhead-2011-manifold-mcmc-published`-kopian, inte den skadade första kopian.

<!-- FILE_INVENTORY -->

## 11. Fullständig filförteckning

Extraktion: Poppler pdftotext 26.04.0, `-layout -enc UTF-8`. Sidgränser bevaras som form feed (`\f`). OCR: Tesseract 5.5.2, engelska, bilder renderade med PyMuPDF i skala 2.5. Originalextraktionerna har inte ersatts eller normaliserats. Hämtad 2026-11-10.

### murray-adams-mackay-2010-elliptical-slice

[PDF](papers/murray-adams-mackay-2010-elliptical-slice.pdf) · [Sökbar råtext](source-text/murray-adams-mackay-2010-elliptical-slice.txt) · [Hämtkälla](https://proceedings.mlr.press/v9/murray10a/murray10a.pdf)

Originalnamn/URL-slut: `murray10a.pdf`. 8 textsegment för 8 PDF-sidor; inga tomma segment. Utvalda relevanta sidor kontrollerade.

PDF SHA-256: `ad6fbbfb8fb7b5899161ab425f0516e81c9c0f3a95ddcd9da04e4745ce232eee`  
Text SHA-256: `8946fb01735bf7717ca594a979bb57cbe3dd8d612374a6af178cbdf5ce2e9654`

### knorr-held-rue-2001-block-updating-preprint

[PDF](papers/knorr-held-rue-2001-block-updating-preprint.pdf) · [Sökbar råtext](source-text/knorr-held-rue-2001-block-updating-preprint.txt) · [Hämtkälla](https://epub.ub.uni-muenchen.de/1600/1/paper_210.pdf)

Originalnamn/URL-slut: `paper_210.pdf`. 39 segment, men råtextens teckenkodning är skadad. Kompletterande OCR, ej ersättning för original.

PDF SHA-256: `8deee3ec205a93abb14de2512844f9fbf3df8eb4fa924bc951af733f22126c7d`  
Text SHA-256: `2c98949e4e408d8ae7c5a1d52d861eefd709b0bf19b49dff1615cb7be32f6585`

[Kompletterande OCR](source-text/knorr-held-rue-2001-block-updating-preprint-ocr.txt), 39 sidor. SHA-256: `d0f31889eda3b0830cfa4272c9904a009e95e0318d9aceb609cf603ee27e002e`. OCR kan förvanska tecken och formler.

### titsias-papaspiliopoulos-2018-auxiliary-gradients-preprint

[PDF](papers/titsias-papaspiliopoulos-2018-auxiliary-gradients-preprint.pdf) · [Sökbar råtext](source-text/titsias-papaspiliopoulos-2018-auxiliary-gradients-preprint.txt) · [Hämtkälla](https://arxiv.org/pdf/1610.09641)

Originalnamn/URL-slut: `1610.09641`. 41 textsegment för 41 PDF-sidor; inga tomma segment. Utvalda relevanta sidor kontrollerade.

PDF SHA-256: `16c1dbf1c295123032b6c576ce680d937055e1717ad58fd41903ef8d8ebfe153`  
Text SHA-256: `1bb23ef283f6d6fb13add474a942a830a0b238a7cb8c37091c9fe70615dc9a9b`

### cotter-etal-2013-mcmc-functions

[PDF](papers/cotter-etal-2013-mcmc-functions.pdf) · [Sökbar råtext](source-text/cotter-etal-2013-mcmc-functions.txt) · [Hämtkälla](https://eprints.maths.manchester.ac.uk/2215/1/1202.0709v3.pdf)

Originalnamn/URL-slut: `1202.0709v3.pdf`. 24 textsegment för 24 PDF-sidor; inga tomma segment. Utvalda relevanta sidor kontrollerade.

PDF SHA-256: `e250115eb80da594ad868dc5799145cc65bbeddaffe7c6a957fd736ba78d768e`  
Text SHA-256: `f02f15468738238e37eeb2fb7059a3f2058b571bcf788e25ad46f43a141537ad`

### beskos-etal-2011-hilbert-hmc-preprint

[PDF](papers/beskos-etal-2011-hilbert-hmc-preprint.pdf) · [Sökbar råtext](source-text/beskos-etal-2011-hilbert-hmc-preprint.txt) · [Hämtkälla](https://www.homepages.ucl.ac.uk/~ucakabe/papers/SPA_11.pdf)

Originalnamn/URL-slut: `SPA_11.pdf`. 35 textsegment för 35 PDF-sidor; inga tomma segment. Utvalda relevanta sidor kontrollerade.

PDF SHA-256: `f83476af99294d660eeea6924a7a201a5bffa8dea4b4c81ec9c2521949004d9f`  
Text SHA-256: `a580d2da9f780eec55e38609e917e24639f6f8ee6178922788f50c464e4b0f14`

### girolami-calderhead-2011-manifold-mcmc

[PDF](papers/girolami-calderhead-2011-manifold-mcmc.pdf) · [Sökbar råtext](source-text/girolami-calderhead-2011-manifold-mcmc.txt) · [Hämtkälla](https://sethna.lassp.cornell.edu/pubPDF/GirolamiComment.pdf)

Originalnamn/URL-slut: `GirolamiComment.pdf`. SKADAD: PDF anger 92 sidor men endast 36 textsegment; tomma segment 12, 18, 26. Används inte för belägg.

PDF SHA-256: `37de680b2a0b7bbed1b976f9cb710d843eedbceacad2f37fdf70787d81427984`  
Text SHA-256: `cdbe7e06f185254fd6ba76545b325ef5f98486692e902519b7503fa1582d3234`

### cui-law-marzouk-2016-dili-preprint

[PDF](papers/cui-law-marzouk-2016-dili-preprint.pdf) · [Sökbar råtext](source-text/cui-law-marzouk-2016-dili-preprint.txt) · [Hämtkälla](https://arxiv.org/pdf/1411.3688)

Originalnamn/URL-slut: `1411.3688`. 42 textsegment för 42 PDF-sidor; inga tomma segment. Utvalda relevanta sidor kontrollerade.

PDF SHA-256: `3de9d4f04d28d090fdbf97b8ed4684705dd7ba7962287a461407db1b63ae20a6`  
Text SHA-256: `184adc4c97b2e00c31450de09ee55157d12a912908a98775dfa99d9b6a45c322`

### papaspiliopoulos-roberts-skold-2007-parametrization

[PDF](papers/papaspiliopoulos-roberts-skold-2007-parametrization.pdf) · [Sökbar råtext](source-text/papaspiliopoulos-roberts-skold-2007-parametrization.txt) · [Hämtkälla](https://arxiv.org/pdf/0708.3797)

Originalnamn/URL-slut: `0708.3797`. 16 textsegment för 16 PDF-sidor; inga tomma segment. Utvalda relevanta sidor kontrollerade.

PDF SHA-256: `21a461334a44211677296e7fc08657ae7381ac80eed1526fb06619abc7c751f1`  
Text SHA-256: `e0247f8b55d233ae0a341c3f1d17bc4577ca0a02189114308ffe0c5cde8caeb5`

### syed-etal-2021-tempering-preprint

[PDF](papers/syed-etal-2021-tempering-preprint.pdf) · [Sökbar råtext](source-text/syed-etal-2021-tempering-preprint.txt) · [Hämtkälla](https://arxiv.org/pdf/1905.02939)

Originalnamn/URL-slut: `1905.02939`. 74 segment. Kraftigt upprepad figurtext; använd selektiv sökning.

PDF SHA-256: `40ddd02b813139413f00a590729d03215021958b0cab041e1c24c893c90083ff`  
Text SHA-256: `e372ae4defd854d57a7fdc5a68680678b58649ad44693034bc27241f6d0a13d6`

### rue-martino-chopin-2008-inla-preprint

[PDF](papers/rue-martino-chopin-2008-inla-preprint.pdf) · [Sökbar råtext](source-text/rue-martino-chopin-2008-inla-preprint.txt) · [Hämtkälla](https://www.statslab.cam.ac.uk/~rjs57/RSS/0708/Rue08.pdf)

Originalnamn/URL-slut: `Rue08.pdf`. 40 textsegment för 40 PDF-sidor; inga tomma segment. Utvalda relevanta sidor kontrollerade.

PDF SHA-256: `0871b0584e8d8bd37f41bbfab9f5e5271c626fa5cad444bf3c447fe712878422`  
Text SHA-256: `ac64c46e9c8a1e5a0b802f151072063d2de230fdfd9669fd18642ea07c24d707`

### andrieu-roberts-2009-pseudo-marginal

[PDF](papers/andrieu-roberts-2009-pseudo-marginal.pdf) · [Sökbar råtext](source-text/andrieu-roberts-2009-pseudo-marginal.txt) · [Hämtkälla](https://arxiv.org/pdf/0903.5480)

Originalnamn/URL-slut: `0903.5480`. 30 textsegment för 30 PDF-sidor; inga tomma segment. Utvalda relevanta sidor kontrollerade.

PDF SHA-256: `4650632179858d85bd745c262ff7c95a26ca9e4cdd25e3ba9c122132c68df484`  
Text SHA-256: `316eda383a881e84ae60767f8a48fce738abecac89160e8108d759b3f667dd71`

### vehtari-etal-2021-rhat-preprint

[PDF](papers/vehtari-etal-2021-rhat-preprint.pdf) · [Sökbar råtext](source-text/vehtari-etal-2021-rhat-preprint.txt) · [Hämtkälla](https://arxiv.org/pdf/1903.08008)

Originalnamn/URL-slut: `1903.08008`. 26 textsegment för 26 PDF-sidor; inga tomma segment. Utvalda relevanta sidor kontrollerade.

PDF SHA-256: `d67a38f7644d7fc77dff8f727c873372f1c43fbb2a507ee146cdc139ae15a24e`  
Text SHA-256: `57d2a2764ad0d168b96093895549bd0f658a5f27a014f4e02c4debfc5fb24e09`

### vehtari-etal-2024-psis

[PDF](papers/vehtari-etal-2024-psis.pdf) · [Sökbar råtext](source-text/vehtari-etal-2024-psis.txt) · [Hämtkälla](https://www.jmlr.org/papers/volume25/19-556/19-556.pdf)

Originalnamn/URL-slut: `19-556.pdf`. 58 textsegment för 58 PDF-sidor; inga tomma segment. Utvalda relevanta sidor kontrollerade.

PDF SHA-256: `2a994fa707dadbf883e238fd14aa93daae2963fedb626c2c656a9cb60fe8e6ab`  
Text SHA-256: `9e8a998b331a87a0a8a218e7ad869ac31f9889c72688e7d559f92241d4a208e0`

### vehtari-gelman-gabry-2017-loo-preprint

[PDF](papers/vehtari-gelman-gabry-2017-loo-preprint.pdf) · [Sökbar råtext](source-text/vehtari-gelman-gabry-2017-loo-preprint.txt) · [Hämtkälla](https://arxiv.org/pdf/1507.04544)

Originalnamn/URL-slut: `1507.04544`. 28 textsegment för 28 PDF-sidor; inga tomma segment. Utvalda relevanta sidor kontrollerade.

PDF SHA-256: `fcb6d1659ab6e2f919372d630f2a8e10347f8737991481cff2791e2fbccb5d97`  
Text SHA-256: `706034c2b9c3d8df0f68ee1a760fe67f4816c8e1ef01acba7d586982725283fc`

### girolami-calderhead-2011-manifold-mcmc-published

[PDF](papers/girolami-calderhead-2011-manifold-mcmc-published.pdf) · [Sökbar råtext](source-text/girolami-calderhead-2011-manifold-mcmc-published.txt) · [Hämtkälla](https://cornebise.com/julien/assets/pdf/girolami2011discussion.pdf)

Originalnamn/URL-slut: `girolami2011discussion.pdf`. 92 textsegment för 92 PDF-sidor; inga tomma segment. Utvalda relevanta sidor kontrollerade.

PDF SHA-256: `87cc93b203eedb51d8a187b0eaf5177495af70983cd5c05746c13f88ebe0403d`  
Text SHA-256: `28ee0b29d3900bc2b00c7e7b6a6af19b87896e970700d69f40a331cd0c44f31b`

