# Mixture LME models with changepoints: Software (R)

<h1>Software to Implement Mixture LME Models with changepoints</h1>

<div class=MsoNormal align=center style='text-align:center'>

<hr size=2 width="100%" align=center>

</div>

<p class=MsoNormal>This directory contains <a href="http://www.r-project.org"><tt><span
style='font-size:10.0pt'>R</span></tt></a> software to implement the
statistical models described in <i><span style='color:black;layout-grid-mode:
line'>Applied Statistic</span></i><span style='color:black;layout-grid-mode:
line'>s, Vol. 53 (3):</span> <span class=spelle><i>Modelling</i></span><i> Growth
and Decline in Lung Function in Duchenne's Muscular Dystrophy with an Augmented
Linear Mixed-Effects Model, </i>by Marc Scott, Robert Norman and Kenneth
Berger.</p>

<p>This research was supported by the National Science Foundation under grant
SES-0088061.<span style='mso-spacerun:yes'>&nbsp; </span>Note: any opinions,
findings, and conclusions or recommendations expressed in this material are
those of the author(s) and do not necessarily reflect the views of the National
Science Foundation</p>

<p>The software can be used freely for non-commercial purposes. You may modify
and distribute the code for non-commercial purposes, as long as
this <a
href="Copyright.html">statement</a> and the contact information is included. </p>

<p>While this software has now been developed and used by the authors, a new
user may well experience some problems or bugs. Please report these by email to
<a href="mailto:marc.scott@nyu.edu">marc.scott@nyu.edu</a>, and we will work
with you to resolve the problem. </p>

<p>The software: </p>

<ul type=circle>
 <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;
     mso-list:l0 level1 lfo1;tab-stops:list .5in'>The <tt><span
     style='font-size:10.0pt'>R</span></tt> <tt><span style='font-size:10.0pt'>code
     </span></tt>consists of a single source file downloadable.&nbsp;&nbsp; <a
     href="http://homepages.nyu.edu/~ms184/Software/mixlme/lmeChgFnsV1.r">Right-click
     here and save code</a>.&nbsp; It should be <tt><span style='font-size:
     10.0pt'>sourced </span></tt>from within <tt><span style='font-size:10.0pt'>R:</span></tt></li>
</ul>

<p><tt><span style='font-size:10.0pt'>source(&quot;lmeChgFnsV1.r&quot;)</span></tt>
</p>

<p>at the prompt. </p>

<p>Minimal documentation is currently available (see below for a rough
guide).&nbsp; The following example code will guide you through a fit of these
models: </p>

<ul type=circle>
 <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;
     mso-list:l2 level1 lfo2;tab-stops:list .5in'>&nbsp;A simulated dataset
     (&quot;d&quot;) is provided for use with the example runs.&nbsp;&nbsp; <a
     href="http://homepages.nyu.edu/~ms184/Software/mixlme/simData.dump">Right-click
     here and save code</a>. (be sure that the file type is .dump, not .txt --
     this can happen in MS Windows environments -- just be sure to change it
     back)&nbsp; It is in &quot;dump&quot; format, so it is <tt><span
     style='font-size:10.0pt'>sourced </span></tt>from within <tt><span
     style='font-size:10.0pt'>R.</span></tt></li>
 <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;
     mso-list:l2 level1 lfo2;tab-stops:list .5in'>&nbsp;This R code loads the
     simulated data, plots it, and then fits one standard lme and one
     changepoint mixture LME model using the data.&nbsp;&nbsp; <a
     href="http://homepages.nyu.edu/~ms184/Software/mixlme/lmechgRuns.r">Right-click
     here and save code for the run</a>. (the plot and summary calls within the
     code work best if called separately using cut &amp; paste, but are
     provided for reference)&nbsp; The runs should be <tt><span
     style='font-size:10.0pt'>sourced </span></tt>from within <tt><span
     style='font-size:10.0pt'>R:</span></tt></li>
</ul>

<p><tt><span style='font-size:10.0pt'>source(&quot;lmechgRuns.r&quot;)</span></tt>
</p>

<p><b><u>Important note:</u></b> the results from the algorithm contain
replications of the original datasource, so proper inference (for the complete
data case) is given by the summary(), rather than print() methods.&nbsp; The
code is currently being modified to provide standard errors using numerical
methods (now available: 11 Nov 04).</p>

<p>The main function that fits the changepoint model is: </p>

<p><span style='font-family:"Courier New"'>lmeChgPt</span>, taking parameters:</p>

<ul type=square>
 <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;
     mso-list:l3 level1 lfo3;tab-stops:list .5in'><span style='font-family:
     "Courier New"'>fixed</span>: fixed effects formula specification</li>
 <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;
     mso-list:l3 level1 lfo3;tab-stops:list .5in'><span style='font-family:
     "Courier New"'>data</span>: data.frame </li>
 <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;
     mso-list:l3 level1 lfo3;tab-stops:list .5in'><span style='font-family:
     "Courier New"'>random</span>: fixed effects formula specification
     (typically of the form Y~X1+X2|Subject)</li>
 <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;
     mso-list:l3 level1 lfo3;tab-stops:list .5in'><span style='font-family:
     "Courier New"'>brkvar</span>: the variable used to specify the breakpoint
     (Age, in the paper)</li>
 <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;
     mso-list:l3 level1 lfo3;tab-stops:list .5in'><span style='font-family:
     "Courier New"'>bkpts</span>: the discrete set of breakpoints to include in
     the model</li>
 <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;
     mso-list:l3 level1 lfo3;tab-stops:list .5in'><span style='font-family:
     "Courier New"'>initProbs</span>: the initial population-level&nbsp;
     probability for each breakpoint</li>
 <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;
     mso-list:l3 level1 lfo3;tab-stops:list .5in'><span style='font-family:
     "Courier New"'>initParams</span>=<span style='font-family:"Courier New"'>list</span>(<span
     style='font-family:"Courier New"'>fixed</span>=, <span style='font-family:
     "Courier New"'>random</span>=), these are the initialization parameters;
     they should follow the structure of fixef() and VarCorr() respectively
     (see worked example)</li>
 <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;
     mso-list:l3 level1 lfo3;tab-stops:list .5in'><span style='font-family:
     "Courier New"'>initChs</span>: a vector of length equal to the number of
     subjects.&nbsp; Initially, each subject must be assigned a specific
     changepoint to start the fitting algorithm.&nbsp; Typically, the
     &quot;middle&quot; changepoint is assigned to every subject.</li>
 <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;
     mso-list:l3 level1 lfo3;tab-stops:list .5in'><span style='font-family:
     "Courier New"'>nRepl</span>: number of times to replicate the data
     (somewhere between 2 and 10 seems to work.&nbsp; We used 6 with good
     success)</li>
 <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;
     mso-list:l3 level1 lfo3;tab-stops:list .5in'><span style='font-family:
     "Courier New"'>max</span>.<span style='font-family:"Courier New"'>fail</span>:
     number of attmepts to make at each iteration.&nbsp; Should a randomization
     (with replication) fail to provide an improvement in likelihood, which can
     happen with the stochastic nature of this EM algorithm, how many times do
     you want to keep trying before moving on to the next iteration (we
     recommend only trying a few times)</li>
 <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;
     mso-list:l3 level1 lfo3;tab-stops:list .5in'><span style='font-family:
     "Courier New"'>max</span>.<span style='font-family:"Courier New"'>iter</span>:
     maximum number of iterations (assuming failure to converge)</li>
 <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;
     mso-list:l3 level1 lfo3;tab-stops:list .5in'><span style='font-family:
     "Courier New"'>thresh</span>: convergence criterion (we use a change in
     loglikelihood of 1e-3)</li>
 <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;
     mso-list:l3 level1 lfo3;tab-stops:list .5in'><span style='font-family:
     "Courier New"'>opt</span>: whether or not to run an optim() call after EM
     finishes (due to instability of optim function on these models, <b>we
     highly recommend setting this to False - see optimizeSE, below</b>)</li>
 <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;
     mso-list:l3 level1 lfo3;tab-stops:list .5in'><span style='font-family:
     "Courier New"'>verbose</span>: indicates whether you wish to view details
     of the EM step (boolean).</li>
</ul>

<p>The function fits the model and returns an object of class <span
style='font-family:"Courier New"'>lmechg</span>, which can be passed to the
support functions given below.</p>

<p><b>Further note:</b> the following three &quot;internal&quot; variables are
generated by the call to lmeChgPt, are part of the replicated data used by the
routines,&nbsp; and may be referenced&nbsp; in formulas passed to support
functions such as <span style='font-family:"Courier New"'>predict</span>:</p>

<ul type=square>
 <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;
     mso-list:l1 level1 lfo4;tab-stops:list .5in'><span style='font-family:
     "Courier New"'>Lo</span>: indicator that observation precedes the break
     point</li>
 <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;
     mso-list:l1 level1 lfo4;tab-stops:list .5in'><span style='font-family:
     "Courier New"'>Hi</span>:&nbsp; indicator that observation is at or
     follows the break point</li>
 <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;
     mso-list:l1 level1 lfo4;tab-stops:list .5in'><span style='font-family:
     "Courier New"'>brk</span>: the break point associated with this subject in
     this replication.&nbsp; This can be used to enforce continuity at the
     breakpoint--see the paper for more details.</li>
</ul>

<p>The following set of &quot;lme-like&quot; functions have been written for
the lmechg object class:</p>

<ul type=circle>
 <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;
     mso-list:l4 level1 lfo5;tab-stops:list .5in'><span style='font-family:
     "Courier New"'>summary</span> (very important - this is best way to view
     results, and generates confidence intervals for many parameters)</li>
 <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;
     mso-list:l4 level1 lfo5;tab-stops:list .5in'><span style='font-family:
     "Courier New"'>print</span></li>
 <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;
     mso-list:l4 level1 lfo5;tab-stops:list .5in'><span style='font-family:
     "Courier New"'>fixef</span> (returns fixed effects)</li>
 <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;
     mso-list:l4 level1 lfo5;tab-stops:list .5in'><span style='font-family:
     "Courier New"'>ranef</span> (returns random effects)</li>
 <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;
     mso-list:l4 level1 lfo5;tab-stops:list .5in'><span style='font-family:
     "Courier New"'>mixProbs</span>.<span style='font-family:"Courier New"'>lmechg</span>
     (yield mixture probabilities) </li>
 <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;
     mso-list:l4 level1 lfo5;tab-stops:list .5in'><span style='font-family:
     "Courier New"'>predict</span> (takes additional argument: popn, a boolean,
     which defaults to False, specifying that the population level mixture
     probabilities should be used in computing individual-specific fixed or
     random effects predictions (BLUPs).&nbsp; Used by default when predictions
     are requested for a subject that was not in the data used to fit the
     model.)</li>
 <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;
     mso-list:l4 level1 lfo5;tab-stops:list .5in'><span style='font-family:
     "Courier New"'>optimizeSE</span> (takes a changepoint model object and
     returns and optimize object not meant for direct use; see <span
     style='font-family:"Courier New"'>updateFixedSE</span> function; this is a
     way to run the optim() call separately from the lmeChgPt call, avoiding
     crashes due to optim() instability)</li>
 <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;
     mso-list:l4 level1 lfo5;tab-stops:list .5in'><span style='font-family:
     "Courier New"'>updateFixedSE</span> (takes a changepoint model object and
     an optimize object and merges these to produce a new changepoint model so
     that summary(new_object) uses results from optim run in calculating
     standard errors of fixed effects.&nbsp; </li>
</ul>

<div class=MsoNormal align=center style='text-align:center'>

<hr size=2 width="100%" align=center>

</div>

<p>Description of these models and an application to growth and decline in lung
function can be found in the paper. </p>

</div>

</body>

</html>
