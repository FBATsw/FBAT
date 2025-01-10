<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>FBAT Team</title>
<meta name="viewport" content="width=device-width, initial-scale=1.0" />
<meta name="description" content="FBAT Team Information" />
<meta name="author" content="FBAT Team" />
<!-- css -->
<link href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" rel="stylesheet">
<link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css" rel="stylesheet">
<style>
  body {
    font-family: Arial, sans-serif;
    background-color: #f9f9f9;
    margin: 0;
    padding: 0;
  }
  #header {
    background-color: #2c3e50;
    color: white;
    padding: 30px 0;
    text-align: center;
  }
  #header h1 {
    margin: 0;
    font-size: 2em;
  }
  #header p {
    margin: 0;
    font-size: 1.2em;
  }
  #content {
    padding: 20px;
    text-align: center;
  }
  #content h2 {
    font-size: 1.8em;
    margin-top: 0;
  }
  #content p, #content ul {
    font-size: 1em;
    line-height: 1.6em;
    text-align: left;
    max-width: 800px;
    margin: 0 auto;
  }
  footer {
    background-color: #2c3e50;
    color: white;
    padding: 10px 0;
    text-align: center;
    font-size: 0.9em;
  }
</style>
</head>
<body>
<div id="header">
  <h1>FBAT TEAM</h1>
  <p>(alphabetical order)<br>
  Gourab De, Julian Hecker, Steve Horvath, Nan Laird, Steve Lake, Christoph Lange, Sharon Lutz, Kristel Van Steen, Adam Sykes, Lin Wang, Wai-Ki Yip, Xin Xu, and Jin Jin Zhou - FBAT</p>
</div>

<div id="content">
  <h2>FBAT</h2>
  <p>FBAT provides the software implementation for several family-based association tests. The FBAT package was developed in the Department of Biostatistics at the Harvard T.H. Chan School of Public Health.</p>
  <ul>
    <li><strong>FBAT v208</strong> – Family Based Association Testing software (Linux executable, 64 bit, region-based extension "gen_rv", development version)</li>
    <li><strong>FBAT v204</strong> – Family Based Association Testing software (Linux executable, 64 bit, before region-based extension, stable version)</li>
    <li><strong>FBAT User's Manual</strong> (current version does not cover new region-based methodology)</li>
    <li><strong>fbat-tour</strong> -- example files I</li>
    <li><strong>FBAT example ped file</strong> -- example file containing 1000 trios and 30 variants</li>
    <li><strong>FBAT source code</strong> is available on GitHub: <a href="https://github.com/FBATsw/FBAT" target="_blank">https://github.com/FBATsw/FBAT</a> under the GNU General Public License v3.0.</li>
  </ul>
  <p>FBAT manual needs to be updated to describe the new region-based association methodology. Essentially, the new methodology is implemented in the command "gen_rv". The general structure is the same as for the old Burden test approach "fbat -v0" and "fbat -v1".</p>
  <p><strong>Caution:</strong> The website <a href="https://sites.google.com/view/fbat-web-page" target="_blank">https://sites.google.com/view/fbat-web-page</a> is outdated and will be offline soon.</p>
  <p>Boston, October 2023</p>

  <h2>If you publish results obtained from using FBAT, please cite:</h2>
  <ul>
    <li>Laird NM, Lange C. Family-based designs in the age of large-scale gene-association studies. Nat Rev Genet. 2006 May;7(5):385-94. doi: 10.1038/nrg1839. PMID: 16619052.</li>
    <li>Horvath S, Xu X, Lake SL, Silverman EK, Weiss ST, Laird NM. Family-based tests for associating haplotypes with general phenotype data: application to asthma genetics. Genet Epidemiol. 2004 Jan;26(1):61-9. doi: 10.1002/gepi.10295. PMID: 14691957.</li>
    <li>Hecker J, Xu X, Townes FW, Loehlein Fier H, Corcoran C, Laird N, Lange C. Family-based tests for associating haplotypes with general phenotype data: Improving the FBAT-haplotype algorithm. Genet Epidemiol. 2018 Feb;42(1):123-126. doi: 10.1002/gepi.22094. Epub 2017 Nov 21. PMID: 29159827; PMCID: PMC5774664.</li>
  </ul>
</div>

<footer>
  <p>&copy; 2023 FBAT Team. All rights reserved.</p>
</footer>
</body>
</html>
