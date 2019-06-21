# context("test get pairs by gene for pairscan")
library(here)

test.data.path <- here("tests/testthat/testdata")
file.name <- file.path(test.data.path, "NON_NZO_Reifsnyder_pgm_CAPE_num.csv")
param.file <- file.path(test.data.path, "cape.parameters.yml")
results.path <- file.path("results")
dir.create(results.path, showWarnings = FALSE)

cross <- read.population(file.name)
cross.obj <- cape2mpp(cross)
geno.obj <- cross.obj$geno.obj$geno

data.obj <- Cape$new(
  parameter_file = param.file,
  results_path = results.path,
  pheno = cross.obj$data.obj$pheno,
  chromosome = cross.obj$data.obj$chromosome,
  marker_num = cross.obj$data.obj$marker_num,
  marker_location = cross.obj$data.obj$marker_location,
  geno_names = dimnames(geno.obj),
  geno = geno.obj
)

# the gene.list below is from https://geneweaver.org/viewgenesetdetails/167950
gene.list <- c("Abl2", "Adss", "Apoa2", "Serpinc1", "Atp1a2", "Atp1b1", "Capn2", "Cd247", "Cd48", "F5", "Crp", "D1Pas1", "Dcaf8", "Ephx1", "Fcer1a", "Fcer1g", "Fcgr2b", "Fcgr3", "Fh1", "Gas5", "Glul", "Hlx", "Ifi203", "Ifi204", "Ly9", "Ncf2", "Pbx1", "Prrx1", "Prox1", "Ptgs2", "Eprs", "Xpr1", "Pdc", "Rxrg", "Apcs", "Sele", "Sell", "Selp", "Spna1", "Nhlh1", "Tgfb2", "Fasl", "Tnr", "Usf1", "Lamc2", "Lamc1", "Lamb3", "Lhx4", "Ifi205", "Pou2f1", "Chml", "Ptpn14", "Mpz", "Hsd11b1", "Tnfsf4", "Xcl1", "Soat1", "Pea15a", "Ppox", "Cacna1e", "Olfr16", "Lefty1", "Rgl1", "Traf5", "Kifap3", "Tfb2m", "Kcnj9", "Dhx9", "Enah", "Rgs16", "Rgs8", "Rgs4", "Itpkb", "Psen2", "Cd244", "Nek2", "Kcnk2", "Atf3", "Prdx6", "H3f3a", "Darc", "Degs1", "Rnasel", "Rgs5", "Astn1", "Fmo3", "Rnf2", "Kcnj10", "Pla2g4a", "Mr1", "Angel2", "Ildr2", "Tada1", "Blzf1", "Myoc", "Fbxo28", "Cacybp", "Pfdn2", "Pycr2", "Casq1", "Fmo1", "Tagln2", "Cenpf", "G0s2", "F11r", "Hsd17b7", "Qsox1", "Itln1", "Dedd", "Pex19", "Uap1", "Copa", "Cd84", "Ier5", "Opn3", "Bpnt1", "Adamts4", "Parp1", "Ush2a", "Dnm3", "Kcnh1", "Uhmk1", "Grem2", "Creg1", "Akt3", "Ddr2", "Slc30a1", "Rgs7", "Nr1i3", "Rfwd2", "Esrrg", "Ifi202b", "Sh2d1b1", "Exo1", "Nit1", "Srp9", "Slamf1", "Mixl1", "Atp1a4", "Rabgap1l", "Zfp238", "Slamf6", "Tor3a", "Usp21", "Tlr5", "Hnrnpu", "Vamp4", "Irf6", "Fmn2", "Aldh9a1", "Lmx1a", "Dusp12", "Tbx19", "Prg4", "Ncstn", "Gpa33", "Alyref2", "Sdhc", "Mosc1", "Ufc1", "Nenf", "1810030J14Rik", "Fam36a", "Pcp4l1", "Mgst3", "Prrc2c", "Tsen15", "Swt1", "Tmem206", "Edem3", "Nuf2", "Rrp15", "Mosc2", "Pigc", "3110045C21Rik", "Adck3", "Nvl", "Pigm", "4930523C07Rik", "Gpatch2", "Arpc5", "Ahctf1", "Tiprl", "Efcab2", "Dusp23", "Mpzl1", "Smyd2", "Rab3gap2", "Klhdc9", "1190005F20Rik", "Tatdn3", "1700025G04Rik", "1700016C15Rik", "Fmo2", "Smyd3", "2810025M15Rik", "2810422O20Rik", "Sft2d2", "Cenpl", "Mpc2", "Nos1ap", "Cep170", "Apobec4", "Mettl13", "Pogk", "Pvrl4", "Shcbp1l", "Ccdc19", "1600012P17Rik", "Iars2", "Acbd6", "Lin9", "Angptl1", "Cnih3", "3110040M04Rik", "1700056E22Rik", "1700057K13Rik", "Ankrd45", "Tmco1", "Rd3", "Cep350", "Npl", "Dcaf6", "Scyl3", "1700034H15Rik", "Spata17", "Slamf8", "Tpr", "4930455F23Rik", "4930527J03Rik", "4930562F07Rik", "Slamf7", "4930558K02Rik", "1700009P17Rik", "1700022P22Rik", "Wdr64", "5830403L16Rik", "Slamf9", "Wdr26", "Sdccag8", "Dtl", "Lrrc52", "1700015E13Rik", "Ints7", "Sccpdh", "9430070O13Rik", "Batf3", "Ralgps2", "4921528O07Rik", "1700084C01Rik", "Cnih4", "D730003I15Rik", "Desi2", "Atf6", "Stx6", "Dusp10", "Rgs18", "Mrps14", "Dpt", "Slc19a2", "B4galt3", "Uck2", "Vangl2", "Igsf9", "Fam129a", "Cadm3", "Diexf", "Kmo", "Glt25d2", "Pyhin1", "Gorab", "Lbr", "AI607873", "Trp53bp2", "Suco", "Susd4", "Mael", "Fcrla", "6330403A02Rik", "Sec16b", "Ivns1abp", "Igsf8", "Nphs2", "Fcgr4", "Acbd3", "Capn8", "Sde2", "Tmem63a", "BC003331", "Ndufs2", "Lyplal1", "Camk1g", "Ppp2r5a", "Fmo4", "Traf3ip3", "A130010J15Rik", "Cdc42bpa", "Rcor3", "Pld5", "Dars2", "Pydc3", "Fam5c", "Fam78b", "7530420F21Rik", "Olfml2b", "Fam5b", "Rps6kc1", "Sertad4", "Lefty2", "Rasal2", "Fam20b", "Nmnat2", "Zbtb37", "Syt14", "Dnahc14", "Kctd3", "Hhat", "Vash2", "C130074G19Rik", "Klhl20", "Mfsd7b", "A730013G03Rik", "Cnst", "Lpgat1", "Kif26b", "BC026585", "Nme7", "Adcy10", "Mark1", "Tnn", "Teddm1", "BC034090", "Tnfsf18", "Rcsd1", "Mir194-1", "Mir205", "Mir214", "Mir215", "Fmo6", "Smg7", "Arhgap30", "Tdrd5", "Hmcn1", "Zfp648", "Mettl11b", "Gpr161", "Dusp27", "Slc30a10", "Rc3h1", "Ccdc121", "Gm821", "Nsl1", "5330438I03Rik", "Gm1305", "Aim2", "Olfr218", "Olfr220", "Olfr231", "Olfr414", "Olfr417", "Olfr419", "Olfr420", "Olfr421-ps1", "Olfr424", "Olfr427", "Olfr429", "Olfr430", "Olfr432", "Olfr433", "Olfr1404", "Olfr1406", "Olfr1408", "A430110L20Rik", "Mnda", "Pappa2", "Dnm3os", "Fcrlb", "Tor1aip1", "Tor1aip2", "BC094916", "Fam71a", "9130409I23Rik", "Tomm40l", "BC055324", "Fmo9", "B020018G12Rik", "4922505E12Rik", "Fcrl6", "Mir199a-2", "Fam163a", "Mir350", "Sh2d1b2", "Mir488", "4930500M09Rik", "Gm10530", "Gm10517", "Gm10521", "Gm9931", "Gm9982", "Gm9765", "7420461P10Rik", "Vsig8", "Gpr52", "Gm4846", "Gm4847", "Gm6177", "Rgs21", "Gm16432", "Gm8214", "Gm5531", "Gm7068", "Gm5533", "Vmn1r1", "Gm7897", "Tstd1", "Gm6185", "Gm7694", "Pydc4", "Gm10176", "Gm10518", "Gm11062", "Gm11074", "Gm9530", "Gm2000", "Gm2061", "Mndal", "Gm3837", "Gm3934", "Gm15852", "Scarna3a", "Snora36b", "Snord47", "Mir1927", "Mir1981", "Mir664", "Gm16548", "n-R5s218", "n-R5s219", "n-R5s220", "Gm16608", "Gm16621", "Gm16701", "Gm16721", "Gm17275", "Gm17311", "Gm17323", "Gm17407", "Gm17444", "Gm17520", "Gm17566", "Gm17653", "Gm17698", "Gm17211", "Mir1843b", "Mir3473c", "Mir5117", "Gm19284", "Gm19846", "Gm20305")

data.obj <- select.markers.for.pairscan.by.gene(
  data.obj, geno.obj, gene.list = gene.list,
  bp.buffer = 10000, organism = "mouse"
)

test_that("test that the marker selection method was updated", {
  expect_equal(data.obj$marker_selection_method, "by.gene")
})

# change the reference allele, and check the change
test_that("test that the result set size is correct", {
  expect_lt(20, dim(pairs.which)[1])
  expect_equal(2, dim(pairs.which)[2])
})