(* GENOMIC DATA *)

(* DEFINITIONS *)
(* GERMLINE data contains unnecessary information. after picking the useful columns, this is how they are arranged: *)
ind1Column = 1;
ind2Column = 2;
chrColumn = 3;
IBDStartColumn = 4;
IBDEndColumn = 5;
ind1HapColumn = 6;
ind2HapColumn = 7;
lengthColumn = 8;



(* generate a vector showing how many times each individual is represented in a new bootstrapped sample *)
bootstrapSamplingCount[region_,toSampleAfrAm_,isCollapsed_]:=
Module[{numInRegion,sample,bsVector},
	(* find the total number of individuals in the specified (sub)region *)
	If[TrueQ[isCollapsed],
		numInRegion=Length[regionalDataSortedCollapsed[[region,If[TrueQ[toSampleAfrAm],1,2]]]]
	,
		numInRegion=Length[regionalDataSorted[[region,If[TrueQ[toSampleAfrAm],1,2]]]]
	];

	(* assign each individual a number (1 ... 'numInRegion') and do a random sampling with replacement of the individuals *)
	sample=RandomChoice[Range[numInRegion],numInRegion];

	(* find how many times each individual is sampled *)
	bsVector=Table[Length[Select[sample,#==i&]],{i,numInRegion}];

	Return[bsVector]
]



(* READING FILTERED IBD DATA *)
scale=1.0; (* make sure GERMLINE output is correctly given in cM instead of bp positions *)
threshold=25000; (* this is the value used in filtering of GERMLINE output *)
SetDirectory["/PATH/TO/GERMLINE_OUTPUT/cMcorrected/filtered/"<>ToString[threshold]<>"_0.1_3.0/"];
If[!ValueQ[totalIBDData],
	Print["Reading filtered IBD data from disk...\n"];

	beginTime=TimeUsed[];

	filename=StringJoin["MERGED_filtered_",ToString[threshold],"_3cM_0.1scale_IBD.haps.match"];
	stream=OpenRead[filename];
	data=ReadList[stream,{Word,Word,Word,Word,Number,Number,Number,Word,Word,Number,Number,Word,Number,Number,Number},WordSeparators->{" ","\t"}];
	Close[stream];
	totalIBDData=Map[{StringSplit[#[[2]],"."][[1]],StringSplit[#[[4]],"."][[1]],#[[5]],#[[6]],#[[7]],StringSplit[#[[2]],"."][[2]],StringSplit[#[[4]],"."][[2]],#[[11]]}&,data];
	Remove[data];
	(*If[Length[totalIBDData]\[Equal]0,Print["No IBD data found. Aborting!"];Abort[]]*)

	Print["Total time spent (sec): ",TimeUsed[]-beginTime,"\n"];
,
	Print["Using the IBD data that is already loaded from disk...\n"]
];

(* generating a list of people who are IBD with others based on GERMLINE data *)
allRelatedSamples=totalIBDData[[All,ind1Column]]\[Union]totalIBDData[[All,ind2Column]];



(* SORTING SAMPLES BY REGION AND BY STATE AND ACCORDING TO AFRICAN/EUROPEAN IDENTIFIER *)
Print["Reading and processing regional information for the individuals..."]

regionCodes=
{
	{"1","Northeast","New England Division (ME,NH,VT,MA,RI,CT)","New England"},
	{"2","Northeast","Middle Atlantic Division (NY,NJ,PA)","Middle Atlantic"},
	{"3","Midwest","East North Central Division (OH,IN,IL,MI,WI)","East North Central"},
	{"4","Midwest","West North Central Division (MN,IA,MO,ND,SD,NE,KS)","West North Central"},
	{"5","South","South Atlantic Division (DE,MD,DC,VA,WV,NC,SC,GA,FL)","South Atlantic"},
	{"6","South","East South Central Division (KY,TN,AL,MS)","East South Central"},
	{"7","South","West South Central Division (AR,LA,OK,TX)","West South Central"},
	{"8","West","Mountain Division (MT,ID,WY,CO,NM,AZ,UT,NV)","Mountain"},
	{"9","West","Pacific Division (WA,OR,CA,AK,HI)","Pacific"},
	{"10","US","N/A"},
	{"11","Foreign Country;Puerto Rico;Other US Territory","N/A"},
	{"98","DK","N/A"},
	{"99","N/A","N/A"},
	{"N/A","N/A","Never interviewed"}
};
StatesByRegions=
{
	{"Connecticut","Maine","Massachusetts","New Hampshire","Rhode Island","Vermont"},
	{"New Jersey","New York","Pennsylvania"},
	{"Illinois","Indiana","Michigan","Ohio","Wisconsin"},
	{"Iowa","Kansas","Minnesota","Missouri","Nebraska","North Dakota","South Dakota"},
	{"Delaware","District Of Columbia","Florida","Georgia","Maryland","North Carolina","South Carolina","Virginia","West Virginia"},
	{"Alabama","Kentucky","Mississippi","Tennessee"},
	{"Arkansas","Louisiana","Oklahoma","Texas"},
	{"Arizona","Colorado","Idaho","Montana","Nevada","New Mexico","Utah","Wyoming"},
	{(*"Alaska",*)"California",(*"Hawaii",*)"Oregon","Washington"}
};



databaseHRSSubjectInfo="/PATH/TO/HRS_METADATA/allIndivs_filtered.BirthSchool2010Region.txt";
HRS\[LetterSpace]Info=Import[databaseHRSSubjectInfo,"CSV"];

databaseSCCSSubjectState="/PATH/TO/SCCS_METADATA/ZipcodeState.txt";
SCCSstateData=Drop[StringSplit[Import[databaseSCCSSubjectState,"List"], "\t"],1];
SCCSregionalData=Map[{#[[1]],Flatten[Position[StatesByRegions,#[[2]]]][[1]]}&,SCCSstateData];

databaseASWSubject="/PATH/TO/1000GENOMES_ASW_METADATA/refpanel.list";
ASWstateData={#[[1]],"Oklahoma"}&/@StringSplit[Import[databaseASWSubject,"List"]]; (* we have assigned them all to Oklahoma, i.e., to West South Central Division *)
ASWregionalData=Map[{#[[1]],Flatten[Position[StatesByRegions,#[[2]]]][[1]]}&,ASWstateData];



(* of the region codes defined above, we focus on the well-defined ones (i.e., 1-9 only) *)
(* the variable subjectsInRegion contains the subjects grouped based on their region of birth *)
(* NOTE: we don't know yet (as of July 30, 2014) where exactly the ASW's in the 1000 Genomes Projects were sequenced at. so their IBD information will not be included in the analysis.
UPDATE: as of September 28, 2014, we will assign ASWs to Oklahoma, i.e., to West South Central Division (see above) *)
subjectsInRegion={};

collapsed=False;
For[region=1,region<=9,region++,(* sifting through HRS sample data *)
	AppendTo[subjectsInRegion,(ToString[#]&/@(Select[HRS\[LetterSpace]Info,#[[(*7*)(*region of birth*)(*9*)(*region of school/10-year-old*)11(*2010 region*)]]==region&][[All,1]](* this include ALL races (white, black, other) *)))\[Intersection]allRelatedSamples]
]

For[region=1,region<=9,region++,(* adding SCCS sample data *)
	SCCSsamplesInRegion=(Select[SCCSregionalData,#[[2]]==region&][[All,1]])\[Intersection]allRelatedSamples;
	If[Length[SCCSsamplesInRegion]!=0,subjectsInRegion[[region]]=Join[subjectsInRegion[[region]],SCCSsamplesInRegion]]
]

For[region=1,region<=9,region++,(* adding ASW sample data *)
	ASWsamplesInRegion=(Select[ASWregionalData,#[[2]]==region&][[All,1]])\[Intersection]allRelatedSamples;
	If[Length[ASWsamplesInRegion]!=0,subjectsInRegion[[region]]=Join[subjectsInRegion[[region]],ASWsamplesInRegion]]
]

collapsed=True; (* collapsing neighboring regions with low number of individuals, (1+2) and (8+9) *)
subjectsInRegionCollapsed={{},Join[subjectsInRegion[[1]],subjectsInRegion[[2]]],subjectsInRegion[[3]],subjectsInRegion[[4]],subjectsInRegion[[5]],subjectsInRegion[[6]],subjectsInRegion[[7]],Join[subjectsInRegion[[8]],subjectsInRegion[[9]]],{}};

(* the variable isAfrAm is a tag (True/False) for the subjects in each region indicating whether they have self-identified as African-Americans or not *)
(* we are using only non-Hispanic African-Americans who are identified correctly in both dbGaP and HRSTracker and were born and currently live within mainland US *)
subjectsAfrAm=ToString[#]&/@(Select[HRS\[LetterSpace]Info,(#[[4]]!="Hispanic"&&#[[6]]=="Black"&&#[[15]]=="AfrAm"&&1<=#[[7]]<=9&&1<=#[[11]]<=9)&][[All,1]]); (* from HRS *)
subjectsAfrAm=Join[subjectsAfrAm,SCCSregionalData[[All,1]]]; (* from SCCS *)
subjectsAfrAm=Join[subjectsAfrAm,ASWregionalData[[All,1]]]; (* from ASW *)
isAfrAm={};
isAfrAmCollapsed={};
For[region=1,region<=9,region++,
	AppendTo[isAfrAm,MemberQ[subjectsAfrAm,#]&/@subjectsInRegion[[region]]];
	AppendTo[isAfrAmCollapsed,MemberQ[subjectsAfrAm,#]&/@subjectsInRegionCollapsed[[region]]]
]

(* we are using only non-Hispanic European-Americans who are identified correctly in both dbGaP and HRSTracker and were born and currently live within mainland US *)
subjectsEurAm=ToString[#]&/@(Select[HRS\[LetterSpace]Info,(#[[4]]!="Hispanic"&&#[[6]]=="White"&&#[[15]]=="Not_AfrAm"&&1<=#[[7]]<=9&&1<=#[[11]]<=9)&][[All,1]]); (* from HRS; no European-Americans in SCCS and ASW *)
isEurAm={};
isEurAmCollapsed={};
For[region=1,region<=9,region++,
	AppendTo[isEurAm,MemberQ[subjectsEurAm,#]&/@subjectsInRegion[[region]]];
	AppendTo[isEurAmCollapsed,MemberQ[subjectsEurAm,#]&/@subjectsInRegionCollapsed[[region]]]
]

(* list of African-Americans in each region *)
AfrAmsInRegion=Table[subjectsInRegion[[i,Flatten[Position[isAfrAm[[i]],True]]]],{i,9}];
AfrAmsInRegionCollapsed={{},Join[AfrAmsInRegion[[1]],AfrAmsInRegion[[2]]],AfrAmsInRegion[[3]],AfrAmsInRegion[[4]],AfrAmsInRegion[[5]],AfrAmsInRegion[[6]],AfrAmsInRegion[[7]],Join[AfrAmsInRegion[[8]],AfrAmsInRegion[[9]]],{}};

(* list of Europeans in each region *)
EUsInRegion=Table[subjectsInRegion[[i,Flatten[Position[isEurAm[[i]],True]]]],{i,9}];
EUsInRegionCollapsed={{},Join[EUsInRegion[[1]],EUsInRegion[[2]]],EUsInRegion[[3]],EUsInRegion[[4]],EUsInRegion[[5]],EUsInRegion[[6]],EUsInRegion[[7]],Join[EUsInRegion[[8]],EUsInRegion[[9]]],{}};

(* regions are sorted from 1 ... 9, and in each region, African-Americans are listed first, then the Europeans: {{{AfrAm_1}, {EU_1}}_1, {{AfrAm_2}, {EU_2}}_2, ..., {{AfrAm_9}, {EU_9}}_9} *)
regionalDataSorted=Table[{AfrAmsInRegion[[i]],EUsInRegion[[i]]},{i,9}];
(* structure of the collapsed list is as follows {{}_1, {{AfrAm_1 + AfrAm_2}, {EU_1 + EU_2}}_2, {{AfrAm_3}, {EU_3}}_3, ..., {{AfrAm_7}, {EU_7}}_7, {{AfrAm_8 + AfrAm_9}, {EU_8 + EU_9}}_8, {}_9} *)
regionalDataSortedCollapsed={{},{Join[regionalDataSorted[[1,1]],regionalDataSorted[[2,1]]],Join[regionalDataSorted[[1,2]],regionalDataSorted[[2,2]]]},regionalDataSorted[[3]],regionalDataSorted[[4]],regionalDataSorted[[5]],regionalDataSorted[[6]],regionalDataSorted[[7]],{Join[regionalDataSorted[[8,1]],regionalDataSorted[[9,1]]],Join[regionalDataSorted[[8,2]],regionalDataSorted[[9,2]]]},{}};

subjectIDs=Flatten[regionalDataSorted];
subjectIDsCollapsed=Flatten[regionalDataSortedCollapsed];

Print["Number of samples in each region: ",Table[Length[subjectsInRegion[[i]]],{i,1,9}]]
Print["Total number of samples: ",Total[Table[Length[subjectsInRegion[[i]]],{i,1,9}]],"\n"]
Print["Number of African-Americans in each region: ",Table[Length[Select[isAfrAm[[i]],#==True&]],{i,1,9}]]
Print["Number of African-Americans in each collapsed region: ",Table[Length[Select[isAfrAmCollapsed[[i]],#==True&]],{i,1,9}]]
Print["Total number of African-Americans: ",Total[Table[Length[Select[isAfrAm[[i]],#==True&]],{i,1,9}]],"\n"]
Print["Number of European-Americans in each region: ",Table[Length[Select[isEurAm[[i]],#==True&]],{i,1,9}]]
Print["Number of European-Americans in each collapsed region: ",Table[Length[Select[isEurAmCollapsed[[i]],#==True&]],{i,1,9}]]
Print["Total number of European-Americans: ",Total[Table[Length[Select[isEurAm[[i]],#==True&]],{i,1,9}]],"\n"]

Print["Creating hash tables for individuals...\n"]
isInCollapsedRegions=Dispatch[(#->True)&/@subjectIDsCollapsed]; (* make a dispatch table to store whether individuals are in the collapsed regions or not *)
getIndex=Dispatch[(#->Flatten[Position[subjectIDsCollapsed,#]][[1]])&/@subjectIDsCollapsed]; (* make a dispatch table to store positions of individuals in the sorted matrix *)



(* READING BIN INFORMATION *)
(* filter IBD data (segments larger than a certain cutoff or segments in a certain length range) *)
Remove[bins,num]

fields[Dynamic[n_]]:=Module[{},Dynamic[Column[Row[{InputField[Dynamic[bins[#,1]],Number],"   ",InputField[Dynamic[bins[#,2]],Number]}]&/@Range[n]]]]

(* having the parentheses is necessary to transform the whole block into one single compound statement; otherwise, the line-breaks would cause Goto[] to not find Label[] *)
(
	Label["getLengthBins"];
	
	bins[1,1]=0;bins[1,2]=\[Infinity];
	
	DialogInput[
	{
		Style["To consider all IBD segments in one bin, press OK without changing anything.",13],Style["NOTE: bins cannot overlap and have to be ordered.",Red,13]
		,
		Row[{"Number of bins for IBD length: ",PopupMenu[Dynamic[num],Range[10]]}]
		,
		Row[{"From:                                                   To:"}]
		,
		DynamicModule[{n=Dynamic[num]},fields[n]]
		,
		DefaultButton[]
	}
	,(*Modal\[Rule]True,*)WindowTitle->"Provide binning information for IBD segments..."
	];
	
	For[i=1,i<=num-1,i++,
		If[bins[i,2]>bins[i+1,1],Goto["getLengthBins"]];
		If[bins[i,2]==bins[i+1,1],bins[i,2]-=10^-3]
	]
)

temp=DownValues[bins];
intervals={};
For[i=1,i<=num,i++,
	AppendTo[intervals,{N[temp[[2 i-1,2]]],N[temp[[2 i,2]]]}]
]
Remove[temp,bins]

intervals



(* PROCESSING AND BINNING IBD SEGMENTS AND CALCULATING RELATEDNESS MATRICES  *)
Print["Creating an array containing IBD individuals and the total length of IBD segments shared..."]

ClearAll[IBDLengthArray,IBDNumberArray]

beginTime=TimeUsed[];
timeElapsed=TimeUsed[];

count=0;
Scan[
(
	ind1=#[[ind1Column]];
	ind2=#[[ind2Column]];
	IBD=(#[[IBDEndColumn]]-#[[IBDStartColumn]])/scale;
	
	isInBins=(#[[1]]<=IBD<=#[[2]])&/@intervals;
	If[MemberQ[isInBins,True],
		binNumber=Flatten[Position[isInBins,True]][[1]];
		
		If[(TrueQ[ind1/.isInCollapsedRegions]&&TrueQ[ind2/.isInCollapsedRegions])&&(ind1!=ind2),
			ind1Index=ind1/.getIndex;
			ind2Index=ind2/.getIndex;
			posInMatrix={Min[ind1Index,ind2Index],Max[ind1Index,ind2Index]};
			
			If[NumberQ[IBDLengthArray[binNumber,posInMatrix]],IBDLengthArray[binNumber,posInMatrix]+=IBD,IBDLengthArray[binNumber,posInMatrix]=IBD];
			If[NumberQ[IBDNumberArray[binNumber,posInMatrix]],IBDNumberArray[binNumber,posInMatrix]+=1,IBDNumberArray[binNumber,posInMatrix]=1]
		]
	];
	
	count++;
	If[Mod[count,250000]==0,
		Print[count,", time spent (sec): ",TimeUsed[]-timeElapsed];
		timeElapsed=TimeUsed[]
	]
)&
,totalIBDData]

Print["Total time spent (sec): ",TimeUsed[]-beginTime,"\n"];

Print["Creating a list of rules, based on the above array, to populate the matrix of IBD relatedness of all individuals..."]

ClearAll[IBDLengthRules,IBDNumberRules]

beginTime=TimeUsed[];

IBDLengthRules=
Map[
(
	tempString=StringDrop[ToString[#],StringLength["HoldPattern["]];
	
	beginBin=First[Flatten[StringPosition[tempString,"["]]]+1;
	firstComma=StringPosition[tempString,","][[1,1]];
	beginPos=First[Flatten[StringPosition[tempString,"{"]]]+1;
	secondComma=StringPosition[tempString,","][[2,1]];
	endPos=First[Flatten[StringPosition[tempString,"}"]]]-1;
	
	binNumber=ToExpression[StringTake[tempString,{beginBin,firstComma-1}]];
	firstPosElement=ToExpression[StringTake[tempString,{beginPos,secondComma-1}]];
	secondPosElement=ToExpression[StringTake[tempString,{secondComma+1,endPos}]];
	{binNumber,{firstPosElement,secondPosElement}}->#[[2]]
)&
,DownValues[IBDLengthArray]];

IBDNumberRules=
Map[
(
	tempString=StringDrop[ToString[#],StringLength["HoldPattern["]];
	
	beginBin=First[Flatten[StringPosition[tempString,"["]]]+1;
	firstComma=StringPosition[tempString,","][[1,1]];
	beginPos=First[Flatten[StringPosition[tempString,"{"]]]+1;
	secondComma=StringPosition[tempString,","][[2,1]];
	endPos=First[Flatten[StringPosition[tempString,"}"]]]-1;
	
	binNumber=ToExpression[StringTake[tempString,{beginBin,firstComma-1}]];
	firstPosElement=ToExpression[StringTake[tempString,{beginPos,secondComma-1}]];
	secondPosElement=ToExpression[StringTake[tempString,{secondComma+1,endPos}]];
	{binNumber,{firstPosElement,secondPosElement}}->#[[2]]
)&
,DownValues[IBDNumberArray]];

Print["Total time spent (sec): ",TimeUsed[]-beginTime,"\n"];

Print["Creating the matrices of IBD relatedness between all individuals..."]

beginTime=TimeUsed[];

IBDLengthBins=
Table[
	rulesPerBin=(#[[1,2]]->#[[2]])&/@Select[IBDLengthRules,#[[1,1]]==binNumber&];
	If[!MemberQ[rulesPerBin[[All,1]],{Length[subjectIDsCollapsed],Length[subjectIDsCollapsed]}],AppendTo[rulesPerBin,({Length[subjectIDsCollapsed],Length[subjectIDsCollapsed]}->0.0)]];
	SparseArray[rulesPerBin]
,{binNumber,num}];

IBDNumberBins=
Table[
	rulesPerBin=(#[[1,2]]->#[[2]])&/@Select[IBDNumberRules,#[[1,1]]==binNumber&];
	If[!MemberQ[rulesPerBin[[All,1]],{Length[subjectIDsCollapsed],Length[subjectIDsCollapsed]}],AppendTo[rulesPerBin,({Length[subjectIDsCollapsed],Length[subjectIDsCollapsed]}->0)]];
	SparseArray[rulesPerBin]
,{binNumber,num}];

IBDShadowBins=Unitize[#]&/@IBDLengthBins;

Print["Total time spent (sec): ",TimeUsed[]-beginTime,"\n"];

Table[Export["MERGED_matrix_IBDLength_"<>ToString[threshold]<>"_bin"<>ToString[binNumber]<>".haps.rua",IBDLengthBins[[binNumber]],"HarwellBoeing"],{binNumber,num}]
Table[Export["MERGED_matrix_IBDNumber_"<>ToString[threshold]<>"_bin"<>ToString[binNumber]<>".haps.rua",IBDNumberBins[[binNumber]],"HarwellBoeing"],{binNumber,num}]



(* REMOVING RELATED INDIVIDUALS *)
(* processing related individuals *)
Print["Removing the contributions of related individuals...\n"]

(* from HRS *)
relatedIndivs={ToString[#[[2]]],ToString[#[[3]]]}&/@Drop[Import["/PATH/TO/HRS_DATA/Kinship_coefficient_table.csv"],1];
(* find their positions in the IBD relatedness matrix *)
relatedIndivsIndex={If[MemberQ[subjectIDsCollapsed,#[[1]]],#[[1]]/.getIndex],If[MemberQ[subjectIDsCollapsed,#[[2]]],#[[2]]/.getIndex]}&/@relatedIndivs;
relatedIndivsIndex=Select[relatedIndivsIndex,!(#[[1]]===Null||#[[2]]===Null)&];
(* remove each pair's respective IBD contribution from the 'Length', 'Number', and 'Shadow' matrices *)
Print["Number of related pairs in HRS: ",Length[relatedIndivsIndex]];
For[i=1,i<=Length[relatedIndivsIndex],i++,
	(* NOTE: matrices are upper-triangular now *)
	ind1=Min[relatedIndivsIndex[[i,1]],relatedIndivsIndex[[i,2]]];
	ind2=Max[relatedIndivsIndex[[i,1]],relatedIndivsIndex[[i,2]]];
	
	For[binNumber=1,binNumber<=num,binNumber++,
		IBDLengthBins[[binNumber]][[ind1,ind2]]=0.0;
		IBDNumberBins[[binNumber]][[ind1,ind2]]=0;
		IBDShadowBins[[binNumber]][[ind1,ind2]]=0
	]
]

(* the family structure in HRS has already been taken into account; for SCCS and ASW (plus the inter-cohort structure), we use the fact that the IIDs of these two cohorts start with a letter (whereas those in HRS start with a number). we, thus, remove all kinships between any two IIDs that start with letters. this is not the best way to handle the familial structure, but it's okay for now. *)
(* from SCCS+ASW *)
relatedIndivs=DeleteCases[If[!(NumberQ[ToExpression[StringTake[ToString[#[[1]]],1]]]&&NumberQ[ToExpression[StringTake[ToString[# [[2]]],1]]]),{ToString[#[[1]]],ToString[#[[2]]]}]&/@Import["/PATH/TO/MERGED_DATA/MERGED_relatedIndivs.list"],Null];
relatedIndivsIndex={If[MemberQ[subjectIDsCollapsed,#[[1]]],#[[1]]/.getIndex],If[MemberQ[subjectIDsCollapsed,#[[2]]],#[[2]]/.getIndex]}&/@relatedIndivs;
relatedIndivsIndex=Select[relatedIndivsIndex,!(#[[1]]===Null||#[[2]]===Null)&];
(* remove each pair's respective IBD contribution from the 'Length', 'Number', and 'Shadow' matrices *)
Print["Number of related pairs in SCCS, ASW, and cross-cohort: ",Length[relatedIndivsIndex]];
For[i=1,i<=Length[relatedIndivsIndex],i++,
	(* NOTE: matrices are upper-triangular now *)
	ind1=Min[relatedIndivsIndex[[i,1]],relatedIndivsIndex[[i,2]]];
	ind2=Max[relatedIndivsIndex[[i,1]],relatedIndivsIndex[[i,2]]];
	
	For[binNumber=1,binNumber<=num,binNumber++,
		IBDLengthBins[[binNumber]][[ind1,ind2]]=0.0;
		IBDNumberBins[[binNumber]][[ind1,ind2]]=0;
		IBDShadowBins[[binNumber]][[ind1,ind2]]=0
	]
]

Print["Total length of shared IBD segments between all individuals..."]
Table[Total[IBDLengthBins[[binNumber]],2],{binNumber,num}]
Print["Total number of shared IBD segments between all individuals..."]
Table[Total[IBDNumberBins[[binNumber]],2],{binNumber,num}]



(* HISTOGRAM OF IBD DISTRIBUTION *)
(* separating African-Americans from European-Americans *)
index\[LetterSpace]AfrAms=Flatten[regionalDataSortedCollapsed[[2;;8]][[All,1]]/.getIndex];
index\[LetterSpace]EurAms=Flatten[regionalDataSortedCollapsed[[2;;8]][[All,2]]/.getIndex];
Table[
	temp=IBDNumberBins[[binNumber]]+Transpose[IBDNumberBins[[binNumber]]];
	
	IBD\[LetterSpace]AfrAms=Total[temp[[index\[LetterSpace]AfrAms]],{2}];
	IBD\[LetterSpace]EurAms=Total[temp[[index\[LetterSpace]EurAms]],{2}];
	
	Print[Max[IBD\[LetterSpace]AfrAms]];
	Print[Histogram[IBD\[LetterSpace]AfrAms,Max[IBD\[LetterSpace]AfrAms]+1,ChartStyle->Directive[EdgeForm[None],Blue],Frame->{True,True,False,False},FrameLabel->{"Number of tracts shared with other samples","Number of samples"},BaseStyle->{FontSize->15}]];
	Print[Max[IBD\[LetterSpace]EurAms]];
	Print[Histogram[IBD\[LetterSpace]EurAms,Max[IBD\[LetterSpace]EurAms]+1,ChartStyle->Directive[EdgeForm[None],Red],Frame->{True,True,False,False},FrameLabel->{"Number of tracts shared with other samples","Number of samples"},BaseStyle->{FontSize->15}]];
	
	Print[Histogram[{IBD\[LetterSpace]EurAms,IBD\[LetterSpace]AfrAms},Max[Max[IBD\[LetterSpace]AfrAms]+1,Max[IBD\[LetterSpace]EurAms]+1],ChartStyle->{Directive[EdgeForm[None],Red],Directive[EdgeForm[None],Blue]},ChartLegends->Placed[{Style["European-Americans",15],Style["African-Americans",15]},Above],Frame->{True,True,False,False},FrameLabel->{"Number of tracts shared with other samples","Number of samples"},BaseStyle->{FontSize->15},PlotRangePadding->{None,None}(*,Axes\[Rule]{True,False},ChartLayout\[Rule]"Overlapped"*)]];
	Export["~/github/AfAmpaper/images/IBDhistogram_bin"<>ToString[binNumber]<>".pdf",Histogram[{IBD\[LetterSpace]EurAms,IBD\[LetterSpace]AfrAms},Max[Max[IBD\[LetterSpace]AfrAms]+1,Max[IBD\[LetterSpace]EurAms]+1],ChartStyle->{Directive[EdgeForm[None],Red],Directive[EdgeForm[None],Blue]},ChartLegends->Placed[{Style["European-Americans",15],Style["African-Americans",15]},Above],Frame->{True,True,False,False},FrameLabel->{"Number of tracts shared with other samples","Number of samples"},BaseStyle->{FontSize->15},PlotRangePadding->{None,None}(*,Axes\[Rule]{True,False},ChartLayout\[Rule]"Overlapped"*)]]
,{binNumber,1,num}];
Remove[temp];



(* VISUALIZATION OF THE RELATEDNESS MATRIX *)
Table[MatrixPlot[IBDLengthBins[[binNumber]]],{binNumber,2,num}]
Table[MatrixPlot[IBDNumberBins[[binNumber]]],{binNumber,2,num}]
Table[MatrixPlot[IBDShadowBins[[binNumber]]],{binNumber,2,num}]



(* COLOR-CODING CENSUS REGIONS *)
ClearAll[styles];
StatesByRegions=
{
	{"Connecticut","Maine","Massachusetts","New Hampshire","Rhode Island","Vermont"},
	{"New Jersey","New York","Pennsylvania"},
	{"Illinois","Indiana","Michigan","Ohio","Wisconsin"},
	{"Iowa","Kansas","Minnesota","Missouri","Nebraska","North Dakota","South Dakota"},
	{"Delaware","District Of Columbia","Florida","Georgia","Maryland","North Carolina","South Carolina","Virginia","West Virginia"},
	{"Alabama","Kentucky","Mississippi","Tennessee"},
	{"Arkansas","Louisiana","Oklahoma","Texas"},
	{"Arizona","Colorado","Idaho","Montana","Nevada","New Mexico","Utah","Wyoming"},
	{(*"Alaska",*)"California",(*"Hawaii",*)"Oregon","Washington"}
};
StatesByRegionsEntity=Table[Entity["AdministrativeDivision",{StringReplace[#," "->""],"UnitedStates"}]&/@StatesByRegions[[i]],{i,9}];
styles=
{
	{GeoStyling["OutlineMap",Darker[Blue]],Polygon[Join[StatesByRegionsEntity[[1]],StatesByRegionsEntity[[2]]]]},
	{GeoStyling["OutlineMap",Lighter[Blue,1/3]],Polygon[StatesByRegionsEntity[[3]]]},
	{GeoStyling["OutlineMap",Lighter[Blue,2/3]],Polygon[StatesByRegionsEntity[[4]]]},
	{GeoStyling["OutlineMap",Orange],Polygon[StatesByRegionsEntity[[5]]]},
	{GeoStyling["OutlineMap",Lighter[Orange]],Polygon[StatesByRegionsEntity[[6]]]},
	{GeoStyling["OutlineMap",Lighter[Orange,2/3]],Polygon[StatesByRegionsEntity[[7]]]},
	{GeoStyling["OutlineMap",Lighter[Red]],Polygon[Join[StatesByRegionsEntity[[8]],StatesByRegionsEntity[[9]]]]}
};

(*GeoGraphics[Flatten[styles],GeoBackground\[Rule]None]*)
regionByColor=(GeoGraphics[#,GeoBackground->None,ImageSize->50]&/@styles)[[{7,3,2,1,6,5,4}]]



(* MOVEMENTS BETWEEN REGION OF BIRTH AND 2010 REGION OF RESIDENCE *)
indivsAnalyzed=Select[Drop[HRS\[LetterSpace]Info,1],MemberQ[subjectIDsCollapsed,ToString[#[[1]]]]&];

temp=Table[0,{9},{9}];
For[i=2,i<=8,i++,
	For[j=2,j<=8,j++,
		temp[[i,j]]=Length[Select[indivsAnalyzed,(#[[4]]!="Hispanic"&&#[[6]]=="Black"&&#[[15]]=="AfrAm"&&#[[7]]==i&&#[[11]]==j)&]]
	]
]

temp[[2,2]]+=Length[Select[indivsAnalyzed,(#[[4]]!="Hispanic"&&#[[6]]=="Black"&&#[[15]]=="AfrAm"&&#[[7]]==1&&#[[11]]==1)&]];
temp[[2,2]]+=Length[Select[indivsAnalyzed,(#[[4]]!="Hispanic"&&#[[6]]=="Black"&&#[[15]]=="AfrAm"&&#[[7]]==2&&#[[11]]==1)&]];
For[j=2,j<=8,j++,
	temp[[2,j]]+=Length[Select[indivsAnalyzed,(#[[4]]!="Hispanic"&&#[[6]]=="Black"&&#[[15]]=="AfrAm"&&#[[7]]==1&&#[[11]]==j)&]]
]
temp[[2,8]]+=Length[Select[indivsAnalyzed,(#[[4]]!="Hispanic"&&#[[6]]=="Black"&&#[[15]]=="AfrAm"&&#[[7]]==1&&#[[11]]==9)&]];
temp[[2,8]]+=Length[Select[indivsAnalyzed,(#[[4]]!="Hispanic"&&#[[6]]=="Black"&&#[[15]]=="AfrAm"&&#[[7]]==2&&#[[11]]==9)&]];

temp[[8,8]]+=Length[Select[indivsAnalyzed,(#[[4]]!="Hispanic"&&#[[6]]=="Black"&&#[[15]]=="AfrAm"&&#[[7]]==9&&#[[11]]==9)&]];
temp[[8,8]]+=Length[Select[indivsAnalyzed,(#[[4]]!="Hispanic"&&#[[6]]=="Black"&&#[[15]]=="AfrAm"&&#[[7]]==8&&#[[11]]==9)&]];
For[j=2,j<=8,j++,
	temp[[8,j]]+=Length[Select[indivsAnalyzed,(#[[4]]!="Hispanic"&&#[[6]]=="Black"&&#[[15]]=="AfrAm"&&#[[7]]==9&&#[[11]]==j)&]]
]
temp[[8,2]]+=Length[Select[indivsAnalyzed,(#[[4]]!="Hispanic"&&#[[6]]=="Black"&&#[[15]]=="AfrAm"&&#[[7]]==9&&#[[11]]==1)&]];
temp[[8,2]]+=Length[Select[indivsAnalyzed,(#[[4]]!="Hispanic"&&#[[6]]=="Black"&&#[[15]]=="AfrAm"&&#[[7]]==8&&#[[11]]==1)&]];

temp//MatrixForm
tempPlot=ArrayPlot[Table[0,{7},{7}],Epilog->{MapIndexed[Text[Style[#1,28],Reverse[#2-1/2]]&,Reverse[temp[[{8,4,3,2,7,6,5},{8,4,3,2,7,6,5}]]],{2}]},Mesh->True,Frame->True,FrameStyle->White,FrameTicks->{{None,MapIndexed[{First[#2],#1}&,regionByColor]},{None,MapIndexed[{First[#2],#1}&,regionByColor]}},ImageSize->750]
Export["region-AfrAm-birth_2010.pdf",tempPlot]
Remove[i,j,temp,tempPlot]

temp=Table[0,{9},{9}];
For[i=2,i<=8,i++,
	For[j=2,j<=8,j++,
		temp[[i,j]]=Length[Select[indivsAnalyzed,(#[[4]]!="Hispanic"&&#[[6]]=="White"&&#[[15]]=="Not_AfrAm"&&#[[7]]==i&&#[[11]]==j)&]]
	]
]

temp[[2,2]]+=Length[Select[indivsAnalyzed,(#[[4]]!="Hispanic"&&#[[6]]=="White"&&#[[15]]=="Not_AfrAm"&&#[[7]]==1&&#[[11]]==1)&]];
temp[[2,2]]+=Length[Select[indivsAnalyzed,(#[[4]]!="Hispanic"&&#[[6]]=="White"&&#[[15]]=="Not_AfrAm"&&#[[7]]==2&&#[[11]]==1)&]];
For[j=2,j<=8,j++,
	temp[[2,j]]+=Length[Select[indivsAnalyzed,(#[[4]]!="Hispanic"&&#[[6]]=="White"&&#[[15]]=="Not_AfrAm"&&#[[7]]==1&&#[[11]]==j)&]]
]
temp[[2,8]]+=Length[Select[indivsAnalyzed,(#[[4]]!="Hispanic"&&#[[6]]=="White"&&#[[15]]=="Not_AfrAm"&&#[[7]]==1&&#[[11]]==9)&]];
temp[[2,8]]+=Length[Select[indivsAnalyzed,(#[[4]]!="Hispanic"&&#[[6]]=="White"&&#[[15]]=="Not_AfrAm"&&#[[7]]==2&&#[[11]]==9)&]];

temp[[8,8]]+=Length[Select[indivsAnalyzed,(#[[4]]!="Hispanic"&&#[[6]]=="White"&&#[[15]]=="Not_AfrAm"&&#[[7]]==9&&#[[11]]==9)&]];
temp[[8,8]]+=Length[Select[indivsAnalyzed,(#[[4]]!="Hispanic"&&#[[6]]=="White"&&#[[15]]=="Not_AfrAm"&&#[[7]]==8&&#[[11]]==9)&]];
For[j=2,j<=8,j++,
	temp[[8,j]]+=Length[Select[indivsAnalyzed,(#[[4]]!="Hispanic"&&#[[6]]=="White"&&#[[15]]=="Not_AfrAm"&&#[[7]]==9&&#[[11]]==j)&]]
]
temp[[8,2]]+=Length[Select[indivsAnalyzed,(#[[4]]!="Hispanic"&&#[[6]]=="White"&&#[[15]]=="Not_AfrAm"&&#[[7]]==9&&#[[11]]==1)&]];
temp[[8,2]]+=Length[Select[indivsAnalyzed,(#[[4]]!="Hispanic"&&#[[6]]=="White"&&#[[15]]=="Not_AfrAm"&&#[[7]]==8&&#[[11]]==1)&]];

temp//MatrixForm
tempPlot=ArrayPlot[Table[0,{7},{7}],Epilog->{MapIndexed[Text[Style[#1,28],Reverse[#2-1/2]]&,Reverse[temp[[{8,4,3,2,7,6,5},{8,4,3,2,7,6,5}]]],{2}]},Mesh->True,Frame->True,FrameStyle->White,FrameTicks->{{None,MapIndexed[{First[#2],#1}&,regionByColor]},{None,MapIndexed[{First[#2],#1}&,regionByColor]}},ImageSize->750]
Export["region-EurAm-birth_2010.pdf",tempPlot]
Remove[i,j,temp,tempPlot]



(* ANALYZING MIGRANTS WHO WERE BORN IN REGION 7 AND ARE LIVING IN REGION 4 IN 2010 *)
migrants\[LetterSpace]7to4=Select[indivsAnalyzed,(#[[4]]!="Hispanic"&&#[[6]]=="Black"&&#[[15]]=="AfrAm"&&#[[7]]==7&&#[[11]]==4)&]
Length[migrants\[LetterSpace]7to4]
migrants\[LetterSpace]7to4\[LetterSpace]index=(ToString[#]&/@migrants\[LetterSpace]7to4[[All,1]])/.getIndex

(* start and end index of African-Americans in each region *)
AfrAms\[LetterSpace]index={};
isAfrAmIBDi=True;
For[i=2,i<=8,i++,
	istart=First[regionalDataSortedCollapsed[[i,If[TrueQ[isAfrAmIBDi],1,2]]]]/.getIndex;
	iend=Last[regionalDataSortedCollapsed[[i,If[TrueQ[isAfrAmIBDi],1,2]]]]/.getIndex;
	Print[i,"\t",istart,"\t\t",iend];
	AppendTo[AfrAms\[LetterSpace]index,{istart;;iend}]
]
AfrAms\[LetterSpace]index
Remove[i,istart,iend,isAfrAmIBDi]

(* IBD of the 7-to-4 African-American individuals to other African-Americans in regions 2 to 8 *)
MatrixPlot[#,ColorFunction->"Monochrome",ColorFunctionScaling->False,AspectRatio->1/2]&/@(Transpose[(IBDLengthBins[[3]]+Transpose[IBDLengthBins[[3]]])[[#,migrants\[LetterSpace]7to4\[LetterSpace]index]]]&/@Flatten[AfrAms\[LetterSpace]index])

(* columns are: total IBD length, number of AfrAms in region, average IBD length per pair (one from the 7-to-4 region, one from the region under consideration) *)
MatrixForm[{Total[#,2],Dimensions[#][[2]],Total[#,2]/Dimensions[#][[2]]}&/@(Transpose[(IBDLengthBins[[3]]+Transpose[IBDLengthBins[[3]]])[[#,migrants\[LetterSpace]7to4\[LetterSpace]index]]]&/@Flatten[AfrAms\[LetterSpace]index])]
(* columns are: total IBD tract number, number of AfrAms in region, average IBD tract number per pair (one from the 7-to-4 region, one from the region under consideration) *)
MatrixForm[{Total[#,2],Dimensions[#][[2]],N[Total[#,2]/Dimensions[#][[2]]]}&/@(Transpose[(IBDNumberBins[[3]]+Transpose[IBDNumberBins[[3]]])[[#,migrants\[LetterSpace]7to4\[LetterSpace]index]]]&/@Flatten[AfrAms\[LetterSpace]index])]
(* columns are: total IBD person number, number of AfrAms in region, average IBD persomn number per pair (one from the 7-to-4 region, one from the region under consideration) *)
MatrixForm[{Total[#,2],Dimensions[#][[2]],N[Total[#,2]/Dimensions[#][[2]]]}&/@(Transpose[(IBDShadowBins[[3]]+Transpose[IBDShadowBins[[3]]])[[#,migrants\[LetterSpace]7to4\[LetterSpace]index]]]&/@Flatten[AfrAms\[LetterSpace]index])]
(* it seems, based on the three metrics above, these 21 7-to-4 African-Americans have the most IBD with region 6. so as they moved northward from region 7 to region 4, they took the IBD that region 6 had with region 7 with themselved to region 4. therefore, by considering the 2010 region of residence, region 6 has stronger IBD with region 4 whereas region 7 has weaker IBD now. *)



(* FOR AFRICAN-AMERICANS -- RELATEDNESS BETWEEN REGIONS FOR DIFFERENT BINS -- USING IBD LENGTH *)
Print["Calculating average IBD relatedness w.r.t. the collapsed regions...\n"]

(* generate AfrAm vs AfrAm IBD data *)
Print["...between African-Americans only (for IBD segments of length between "<>ToString[intervals]<>")"]
Print[Style["using average LENGTH",18,Red]]

AfrAmAfrAmIBDLength=Table[Table[0,{9},{9}],{num}];
AfrAmAfrAmIBDNumber=Table[Table[0,{9},{9}],{num}];

isAfrAmIBDi=True;isAfrAmIBDj=True;
For[binNumber=1,binNumber<=num,binNumber++,
	For[i=2,i<=8,i++,
		istart=First[regionalDataSortedCollapsed[[i,If[TrueQ[isAfrAmIBDi],1,2]]]]/.getIndex;
		iend=Last[regionalDataSortedCollapsed[[i,If[TrueQ[isAfrAmIBDi],1,2]]]]/.getIndex;
		
		For[j=i,j<=8,j++,
			jstart=First[regionalDataSortedCollapsed[[j,If[TrueQ[isAfrAmIBDj],1,2]]]]/.getIndex;
			jend=Last[regionalDataSortedCollapsed[[j,If[TrueQ[isAfrAmIBDj],1,2]]]]/.getIndex;
			
			If[i!=j,
				pairs=(iend-istart+1) (jend-jstart+1),
				pairs=(iend-istart+1) (iend-istart)/2
			];
			
			AfrAmAfrAmIBDLength[[binNumber]][[i,j]]=Total[IBDLengthBins[[binNumber,istart;;iend,jstart;;jend]],2]/pairs;
			AfrAmAfrAmIBDNumber[[binNumber]][[i,j]]=Total[IBDShadowBins[[binNumber,istart;;iend,jstart;;jend]],2]
		]
	]
]

Table[MatrixForm[AfrAmAfrAmIBDLength[[binNumber]]],{binNumber,num}]
{Min[Select[Flatten[AfrAmAfrAmIBDLength[[1]]],#!=0&]],Max[AfrAmAfrAmIBDLength[[1]]]}
{Min[Select[Flatten[AfrAmAfrAmIBDLength[[2]]],#!=0&]],Max[AfrAmAfrAmIBDLength[[2]]]}
{Min[Select[Flatten[AfrAmAfrAmIBDLength[[3]]],#!=0&]],Max[AfrAmAfrAmIBDLength[[3]]]}

ArrayPlot[AfrAmAfrAmIBDLength[[#]][[2;;8,2;;8]],Mesh->True,FrameTicks->False,ColorFunctionScaling->False,PlotRange->{Min[Select[Flatten[AfrAmAfrAmIBDLength[[#]]],#!=0&]],Max[AfrAmAfrAmIBDLength[[#]]]}]&/@Range[num]

plots\[LetterSpace]length=
Table[
	temp=AfrAmAfrAmIBDLength[[binNumber]]+Transpose[AfrAmAfrAmIBDLength[[binNumber]]]-DiagonalMatrix[Diagonal[AfrAmAfrAmIBDLength[[binNumber]]]];
	ArrayPlot[Table[0,{7},{7}],Epilog->{MapIndexed[If[#1!=0,Text[Style[NumberForm[#1,{6,3}],28],Reverse[#2-1/2]]]&,Reverse[UpperTriangularize[temp[[{8,4,3,2,7,6,5},{8,4,3,2,7,6,5}]]]],{2}]},Mesh->True,Frame->True,FrameStyle->White,FrameTicks->{{None,MapIndexed[{First[#2],#1}&,regionByColor]},{None,MapIndexed[{First[#2],#1}&,regionByColor]}},ImageSize->750]
,{binNumber,num}]
Table[Export["regionIBD_"<>If[isAfrAmIBDi,"AfrAm","EurAm"]<>"-"<>If[isAfrAmIBDj,"AfrAm","EurAm"]<>"_bin"<>ToString[binNumber]<>"_length.pdf",plots\[LetterSpace]length[[binNumber]]],{binNumber,num}]
Remove[temp]

temp=(AfrAmAfrAmIBDLength[[#]]+Transpose[AfrAmAfrAmIBDLength[[#]]]-DiagonalMatrix[Diagonal[AfrAmAfrAmIBDLength[[#]]]])&/@Range[num];
plots\[LetterSpace]scale=ArrayPlot[UpperTriangularize[temp[[#]][[{8,4,3,2,7,6,5},{8,4,3,2,7,6,5}]]](*temp\[LeftDoubleBracket]#\[RightDoubleBracket]\[LeftDoubleBracket]{8,4,3,2,7,6,5},{8,4,3,2,7,6,5}\[RightDoubleBracket]*),(*Mesh\[Rule]True,*)Frame->True,FrameStyle->White,FrameTicks->{{None,MapIndexed[{First[#2],#1}&,regionByColor]},{None,MapIndexed[{First[#2],#1}&,regionByColor]}},ColorFunctionScaling->False,PlotRange->{Min[Select[Flatten[AfrAmAfrAmIBDLength[[#]]],#!=0&]],Max[AfrAmAfrAmIBDLength[[#]]]},PlotLegends->Placed[BarLegend[Automatic,LabelStyle->{Black,FontSize->30}],Bottom],ImageSize->750]&/@Range[num]
Table[Export["regionIBD_"<>If[isAfrAmIBDi,"AfrAm","EurAm"]<>"-"<>If[isAfrAmIBDj,"AfrAm","EurAm"]<>"_bin"<>ToString[binNumber]<>"_scale.pdf",plots\[LetterSpace]scale[[binNumber]]],{binNumber,num}]
Remove[temp]

plots\[LetterSpace]number=
Table[
	temp=AfrAmAfrAmIBDNumber[[binNumber]]+Transpose[AfrAmAfrAmIBDNumber[[binNumber]]]-DiagonalMatrix[Diagonal[AfrAmAfrAmIBDNumber[[binNumber]]]];
	ArrayPlot[Table[0,{7},{7}],Epilog->{MapIndexed[If[#1!=0,Text[Style[#1,Large],Reverse[#2-1/2]]]&,Reverse[UpperTriangularize[temp[[{8,4,3,2,7,6,5},{8,4,3,2,7,6,5}]]]],{2}]},Mesh->True,Frame->True,FrameStyle->White,FrameTicks->{{None,MapIndexed[{First[#2],#1}&,regionByColor]},{None,MapIndexed[{First[#2],#1}&,regionByColor]}},ImageSize->750]
,{binNumber,num}]
Table[Export["regionIBD_"<>If[isAfrAmIBDi,"AfrAm","EurAm"]<>"-"<>If[isAfrAmIBDj,"AfrAm","EurAm"]<>"_bin"<>ToString[binNumber]<>"_number.pdf",plots\[LetterSpace]number[[binNumber]]],{binNumber,num}]
Remove[temp]



(* FOR EUROPEAN-AMERICANS -- RELATEDNESS BETWEEN REGIONS FOR DIFFERENT BINS -- USING IBD LENGTH *)
(* generate EU vs EU IBD data *)
Print["...between Europeans only (for IBD segments of length between "<>ToString[intervals]<>")"]

EUEUIBDLength=Table[Table[0,{9},{9}],{num}];
EUEUIBDNumber=Table[Table[0,{9},{9}],{num}];

isAfrAmIBDi=False;isAfrAmIBDj=False;
For[binNumber=1,binNumber<=num,binNumber++,
	For[i=2,i<=8,i++,
		istart=First[regionalDataSortedCollapsed[[i,If[TrueQ[isAfrAmIBDi],1,2]]]]/.getIndex;
		iend=Last[regionalDataSortedCollapsed[[i,If[TrueQ[isAfrAmIBDi],1,2]]]]/.getIndex;
		
		For[j=i,j<=8,j++,
			jstart=First[regionalDataSortedCollapsed[[j,If[TrueQ[isAfrAmIBDj],1,2]]]]/.getIndex;
			jend=Last[regionalDataSortedCollapsed[[j,If[TrueQ[isAfrAmIBDj],1,2]]]]/.getIndex;
			
			If[i!=j,
				pairs=(iend-istart+1) (jend-jstart+1),
				pairs=(iend-istart+1) (iend-istart)/2
			];
			
			EUEUIBDLength[[binNumber]][[i,j]]=Total[IBDLengthBins[[binNumber,istart;;iend,jstart;;jend]],2]/pairs;
			EUEUIBDNumber[[binNumber]][[i,j]]=Total[IBDShadowBins[[binNumber,istart;;iend,jstart;;jend]],2]
		]
	]
]

Table[MatrixForm[EUEUIBDLength[[binNumber]]],{binNumber,num}]
{Min[Select[Flatten[EUEUIBDLength[[1]]],#!=0&]],Max[EUEUIBDLength[[1]]]}
{Min[Select[Flatten[EUEUIBDLength[[2]]],#!=0&]],Max[EUEUIBDLength[[2]]]}
{Min[Select[Flatten[EUEUIBDLength[[3]]],#!=0&]],Max[EUEUIBDLength[[3]]]}

ArrayPlot[EUEUIBDLength[[#]][[2;;8,2;;8]],Mesh->True,FrameTicks->False,ColorFunctionScaling->False,PlotRange->{Min[Select[Flatten[EUEUIBDLength[[#]]],#!=0&]],Max[EUEUIBDLength[[#]]]}]&/@Range[num]

plots\[LetterSpace]length=
Table[
	temp=EUEUIBDLength[[binNumber]]+Transpose[EUEUIBDLength[[binNumber]]]-DiagonalMatrix[Diagonal[EUEUIBDLength[[binNumber]]]];
	ArrayPlot[Table[0,{7},{7}],Epilog->{MapIndexed[If[#1!=0,Text[Style[NumberForm[#1,{6,3}],28],Reverse[#2-1/2]]]&,Reverse[UpperTriangularize[temp[[{8,4,3,2,7,6,5},{8,4,3,2,7,6,5}]]]],{2}]},Mesh->True,Frame->True,FrameStyle->White,FrameTicks->{{None,MapIndexed[{First[#2],#1}&,regionByColor]},{None,MapIndexed[{First[#2],#1}&,regionByColor]}},ImageSize->750]
,{binNumber,num}]
Table[Export["regionIBD_"<>If[isAfrAmIBDi,"AfrAm","EurAm"]<>"-"<>If[isAfrAmIBDj,"AfrAm","EurAm"]<>"_bin"<>ToString[binNumber]<>"_length.pdf",plots\[LetterSpace]length[[binNumber]]],{binNumber,num}]
Remove[temp]

temp=(EUEUIBDLength[[#]]+Transpose[EUEUIBDLength[[#]]]-DiagonalMatrix[Diagonal[EUEUIBDLength[[#]]]])&/@Range[num];
plots\[LetterSpace]scale=ArrayPlot[UpperTriangularize[temp[[#]][[{8,4,3,2,7,6,5},{8,4,3,2,7,6,5}]]](*temp\[LeftDoubleBracket]#\[RightDoubleBracket]\[LeftDoubleBracket]{8,4,3,2,7,6,5},{8,4,3,2,7,6,5}\[RightDoubleBracket]*),(*Mesh\[Rule]True,*)Frame->True,FrameStyle->White,FrameTicks->{{None,MapIndexed[{First[#2],#1}&,regionByColor]},{None,MapIndexed[{First[#2],#1}&,regionByColor]}},ColorFunctionScaling->False,PlotRange->{Min[Select[Flatten[EUEUIBDLength[[#]]],#!=0&]],Max[EUEUIBDLength[[#]]]},PlotLegends->Placed[BarLegend[Automatic,LabelStyle->{Black,FontSize->30}],Bottom],ImageSize->750]&/@Range[num]
Table[Export["regionIBD_"<>If[isAfrAmIBDi,"AfrAm","EurAm"]<>"-"<>If[isAfrAmIBDj,"AfrAm","EurAm"]<>"_bin"<>ToString[binNumber]<>"_scale.pdf",plots\[LetterSpace]scale[[binNumber]]],{binNumber,num}]
Remove[temp]

plots\[LetterSpace]number=
Table[
	temp=EUEUIBDNumber[[binNumber]]+Transpose[EUEUIBDNumber[[binNumber]]]-DiagonalMatrix[Diagonal[EUEUIBDNumber[[binNumber]]]];
	ArrayPlot[Table[0,{7},{7}],Epilog->{MapIndexed[If[#1!=0,Text[Style[#1,Large],Reverse[#2-1/2]]]&,Reverse[UpperTriangularize[temp[[{8,4,3,2,7,6,5},{8,4,3,2,7,6,5}]]]],{2}]},Mesh->True,Frame->True,FrameStyle->White,FrameTicks->{{None,MapIndexed[{First[#2],#1}&,regionByColor]},{None,MapIndexed[{First[#2],#1}&,regionByColor]}},ImageSize->750]
,{binNumber,num}]
Table[Export["regionIBD_"<>If[isAfrAmIBDi,"AfrAm","EurAm"]<>"-"<>If[isAfrAmIBDj,"AfrAm","EurAm"]<>"_bin"<>ToString[binNumber]<>"_number.pdf",plots\[LetterSpace]number[[binNumber]]],{binNumber,num}]
Remove[temp]



(* FOR AFRICAN-AMERICANS AND EUROPEAN-AMERICANS -- RELATEDNESS BETWEEN REGIONS FOR DIFFERENT BINS -- USING IBD LENGTH *)
(* generate AfrAm vs EU IBD data *)
Print["...between African-Americans and Europeans (y-axis: AfrAm, x-axis: EurAm) (for IBD segments of length between "<>ToString[intervals]<>")"]
Print[Style["using average LENGTH",18,Red]]

(* the matrices constructed above are upper-triangular (with the diagonal elements being zero); the full IBD matrix is constructed as shown below *) 
l=(#+Transpose[#])&/@IBDLengthBins;
s=(#+Transpose[#])&/@IBDShadowBins; (* NOTE: i was using IBDNumber here before; i believe that was wrong, and IBDShadow should be used instead. *)

AfrAmEUIBDLength=Table[Table[0,{9},{9}],{num}];
AfrAmEUIBDNumber=Table[Table[0,{9},{9}],{num}];

isAfrAmIBDi=True;isAfrAmIBDj=False;
For[binNumber=1,binNumber<=num,binNumber++,
	For[i=2,i<=8,i++,
		istart=First[regionalDataSortedCollapsed[[i,If[TrueQ[isAfrAmIBDi],1,2]]]]/.getIndex;
		iend=Last[regionalDataSortedCollapsed[[i,If[TrueQ[isAfrAmIBDi],1,2]]]]/.getIndex;
		
		For[j=2(*i*),j<=8,j++,(* unlike the other ones above, this matrix is not symmetric! *)
			jstart=First[regionalDataSortedCollapsed[[j,If[TrueQ[isAfrAmIBDj],1,2]]]]/.getIndex;
			jend=Last[regionalDataSortedCollapsed[[j,If[TrueQ[isAfrAmIBDj],1,2]]]]/.getIndex;
			
			pairs=(iend-istart+1) (jend-jstart+1);
			
			AfrAmEUIBDLength[[binNumber]][[i,j]]=Total[l[[binNumber,istart;;iend,jstart;;jend]],2]/pairs;
			AfrAmEUIBDNumber[[binNumber]][[i,j]]=Total[s[[binNumber,istart;;iend,jstart;;jend]],2]
		]
	]
]

Table[MatrixForm[AfrAmEUIBDLength[[binNumber]]],{binNumber,num}]
{Min[Select[Flatten[AfrAmEUIBDLength[[1]]],#!=0&]],Max[AfrAmEUIBDLength[[1]]]}
{Min[Select[Flatten[AfrAmEUIBDLength[[2]]],#!=0&]],Max[AfrAmEUIBDLength[[2]]]}
{Min[Select[Flatten[AfrAmEUIBDLength[[3]]],#!=0&]],Max[AfrAmEUIBDLength[[3]]]}

ArrayPlot[AfrAmEUIBDLength[[#]][[2;;8,2;;8]],Mesh->True,FrameTicks->False,ColorFunctionScaling->False,PlotRange->{Min[Select[Flatten[AfrAmEUIBDLength[[#]]],#!=0&]],Max[AfrAmEUIBDLength[[#]]]}]&/@Range[num]

plots\[LetterSpace]length=
Table[
	ArrayPlot[Table[0,{7},{7}],Epilog->{MapIndexed[If[#1!=0,Text[Style[NumberForm[#1,{6,4}],28],Reverse[#2-1/2]]]&,Reverse[AfrAmEUIBDLength[[binNumber]][[{8,4,3,2,7,6,5},{8,4,3,2,7,6,5}]]],{2}]},Mesh->True,Frame->True,FrameStyle->White,FrameTicks->{{None,MapIndexed[{First[#2],#1}&,regionByColor]},{None,MapIndexed[{First[#2],#1}&,regionByColor]}},ImageSize->750]
,{binNumber,num}]
Table[Export["regionIBD_"<>If[isAfrAmIBDi,"AfrAm","EurAm"]<>"-"<>If[isAfrAmIBDj,"AfrAm","EurAm"]<>"_bin"<>ToString[binNumber]<>"_length.pdf",plots\[LetterSpace]length[[binNumber]]],{binNumber,num}]

plots\[LetterSpace]scale=ArrayPlot[AfrAmEUIBDLength[[#]][[{8,4,3,2,7,6,5},{8,4,3,2,7,6,5}]],(*Mesh\[Rule]True,*)Frame->True,FrameStyle->White,FrameTicks->{{None,MapIndexed[{First[#2],#1}&,regionByColor]},{None,MapIndexed[{First[#2],#1}&,regionByColor]}},ColorFunctionScaling->False,PlotRange->{Min[Select[Flatten[AfrAmEUIBDLength[[#]]],#!=0&]],Max[AfrAmEUIBDLength[[#]]]},PlotLegends->Placed[BarLegend[Automatic,LabelStyle->{Black,FontSize->30}],Bottom],ImageSize->750]&/@Range[num]
Table[Export["regionIBD_"<>If[isAfrAmIBDi,"AfrAm","EurAm"]<>"-"<>If[isAfrAmIBDj,"AfrAm","EurAm"]<>"_bin"<>ToString[binNumber]<>"_scale.pdf",plots\[LetterSpace]scale[[binNumber]]],{binNumber,num}]

plots\[LetterSpace]number=
Table[
	ArrayPlot[Table[0,{7},{7}],Epilog->{MapIndexed[If[#1!=0,Text[Style[#1,Large],Reverse[#2-1/2]]]&,Reverse[AfrAmEUIBDNumber[[binNumber]][[{8,4,3,2,7,6,5},{8,4,3,2,7,6,5}]]],{2}]},Mesh->True,Frame->True,FrameStyle->White,FrameTicks->{{None,MapIndexed[{First[#2],#1}&,regionByColor]},{None,MapIndexed[{First[#2],#1}&,regionByColor]}},ImageSize->750]
,{binNumber,num}]
Table[Export["regionIBD_"<>If[isAfrAmIBDi,"AfrAm","EurAm"]<>"-"<>If[isAfrAmIBDj,"AfrAm","EurAm"]<>"_bin"<>ToString[binNumber]<>"_number.pdf",plots\[LetterSpace]number[[binNumber]]],{binNumber,num}]

(* using a purple pallette (to avoid misinterpretation!) *)
plot3\[LetterSpace]scale\[LetterSpace]PURPLE=ArrayPlot[AfrAmEUIBDLength[[3]][[{8,4,3,2,7,6,5},{8,4,3,2,7,6,5}]],(*Mesh\[Rule]True,*)Frame->True,FrameStyle->White,FrameTicks->{{None,MapIndexed[{First[#2],#1}&,regionByColor]},{None,MapIndexed[{First[#2],#1}&,regionByColor]}},ColorFunction->(Blend[{White,Darker[Purple,1/2]},#]&),(*ColorFunctionScaling\[Rule]False,*)PlotRange->{Min[Select[Flatten[AfrAmEUIBDLength[[3]]],#!=0&]],Max[AfrAmEUIBDLength[[3]]]},PlotLegends->Placed[BarLegend[{Automatic,{10^3 Min[Select[Flatten[AfrAmEUIBDLength[[3]]],#!=0&]],10^3 Max[AfrAmEUIBDLength[[3]]]}}(*,LegendLabel\[Rule]Placed["x10^3",After]*),LabelStyle->{Black,FontSize->50}],Bottom],ImageSize->750]
Export["regionIBD_"<>If[isAfrAmIBDi,"AfrAm","EurAm"]<>"-"<>If[isAfrAmIBDj,"AfrAm","EurAm"]<>"_bin3_scale_PURPLE.pdf",plot3\[LetterSpace]scale\[LetterSpace]PURPLE]



(* MAPS -- FOR AFRICAN-AMERICANS -- RELATEDNESS BETWEEN REGIONS FOR DIFFERENT BINS -- USING IBD LENGTH *)
numAfrAmInRegionCollapsed=Length[#]&/@AfrAmsInRegionCollapsed (* number of African-Americans in each region *)
Total[numAfrAmInRegionCollapsed]

maskAfrAmAfrAm=Table[0,{9},{9}]; (* we will remove region connections if the possible number of pairs of individuals for those two regions fall below 10000 *)
For[i=1,i<=9,i++,
	For[j=1,j<=9,j++,
		If[i!=j,
			numPairsAfrAm=numAfrAmInRegionCollapsed[[i]] numAfrAmInRegionCollapsed[[j]],(* pairs form different regions *)
			numPairsAfrAm=numAfrAmInRegionCollapsed[[i]] (numAfrAmInRegionCollapsed[[i]]-1)/2 (* pairs from the same region *)
		];
		
		If[numPairsAfrAm>=10000,maskAfrAmAfrAm[[i,j]]=1] (* build the mask to apply to the IBD matrix*)
	]
]
maskAfrAmAfrAm//MatrixForm

Print[Style["AfrAm vs AfrAm",18,Blue]]

cutoff=0.25; (* rescaling total IBD length to [cutoff, 1] *)
maxThickness=0.0125;

StatesByRegions=
{
	{"Connecticut","Maine","Massachusetts","New Hampshire","Rhode Island","Vermont"},
	{"New Jersey","New York","Pennsylvania"},
	{"Illinois","Indiana","Michigan","Ohio","Wisconsin"},
	{"Iowa","Kansas","Minnesota","Missouri","Nebraska","North Dakota","South Dakota"},
	{"Delaware","District Of Columbia","Florida","Georgia","Maryland","North Carolina","South Carolina","Virginia","West Virginia"},
	{"Alabama","Kentucky","Mississippi","Tennessee"},
	{"Arkansas","Louisiana","Oklahoma","Texas"},
	{"Arizona","Colorado","Idaho","Montana","Nevada","New Mexico","Utah","Wyoming"},
	{(*"Alaska",*)"California",(*"Hawaii",*)"Oregon","Washington"}
};
StatesByRegionsEntity=Table[Entity["AdministrativeDivision",{StringReplace[#," "->""],"UnitedStates"}]&/@StatesByRegions[[i]],{i,9}];

(*
regionHub=
{
	{"Boston","Massachusetts","UnitedStates"},
	(*{"Philadelphia","Pennsylvania","UnitedStates"}*){"Burlington","Vermont","UnitedStates"},
	{"Urbana","Illinois","UnitedStates"}(*{"Chicago","Illinois","UnitedStates"}*),
	(*{"Omaha","Nebraska","UnitedStates"}*){"Wall","SouthDakota","UnitedStates"}(*{"Fargo","NorthDakota","UnitedStates"}*),
	{"Hinesville","Georgia","UnitedStates"}(*{"Jacksonville","Florida","UnitedStates"}*),
	(*{"Jackson","Mississippi","UnitedStates"}*){"Waynesboro","Mississippi","UnitedStates"},
	(*{"NewOrleans","Louisiana","UnitedStates"}*){"CorpusChristi","Texas","UnitedStates"},
	(*{"Denver","Colorado","UnitedStates"}*)(*{"SaltLakeCity","Utah","UnitedStates"}*){"LasVegas","Nevada","UnitedStates"},
	{"SanFrancisco","California","UnitedStates"}
};
*) (* these are the old hubs -- i prefer the new ones below *)
regionHub=
{
	{"Boston","Massachusetts","UnitedStates"},
	(*{"Philadelphia","Pennsylvania","UnitedStates"}*){"Burlington","Vermont","UnitedStates"},
	{"Chicago","Illinois","UnitedStates"},
	(*{"Omaha","Nebraska","UnitedStates"}*)(*{"Minneapolis","Minnesota","UnitedStates"}*){"Fargo","NorthDakota","UnitedStates"},
	(*{"Atlanta","Georgia","UnitedStates"}*){"Jacksonville","Florida","UnitedStates"},
	(*{"Jackson","Mississippi","UnitedStates"}*){"Waynesboro","Mississippi","UnitedStates"},
	(*{"NewOrleans","Louisiana","UnitedStates"}*){"Houston","Texas","UnitedStates"},
	(*{"Denver","Colorado","UnitedStates"}*){"SaltLakeCity","Utah","UnitedStates"},
	{"SanFrancisco","California","UnitedStates"}
};
regionHubEntity=Entity["City",#]&/@regionHub;

Table[
	Print[Style["\n\n\nBIN "<>ToString[binNumber]<>": "<>ToString[intervals[[binNumber]]]<>"cM",26,Red]];
	
	interRegionIBD=(AfrAmAfrAmIBDLength[[binNumber]]+Transpose[AfrAmAfrAmIBDLength[[binNumber]]]-2 DiagonalMatrix[Diagonal[AfrAmAfrAmIBDLength[[binNumber]]]]) maskAfrAmAfrAm(* apply the mask constructed above *);
	minIBD=Min[Select[Flatten[interRegionIBD],#!=0&]];
	maxIBD=Max[interRegionIBD];
	(* fixing the scaling of the lines for AfrAm-AfrAm and EurAm-EurAm plots to the same scale *)
	(*interRegionIBDEUEU=(EUEUIBDLength\[LeftDoubleBracket]binNumber\[RightDoubleBracket]+Transpose[EUEUIBDLength\[LeftDoubleBracket]binNumber\[RightDoubleBracket]]-2 DiagonalMatrix[Diagonal[EUEUIBDLength\[LeftDoubleBracket]binNumber\[RightDoubleBracket]]]);
	minIBD=Min[Min[Select[Flatten[interRegionIBD],#\[NotEqual]0&]],Min[Select[Flatten[interRegionIBDEUEU],#\[NotEqual]0&]]];
	maxIBD=Max[Max[interRegionIBD],Max[interRegionIBDEUEU]];*)
	normalized=Map[If[#!=0,(1-cutoff) (#-minIBD)/(maxIBD-minIBD)+cutoff,0]&,interRegionIBD,{2}];
	
	ClearAll[styles];
	styles=
	{
		{GeoStyling["OutlineMap",Darker[Blue]],Polygon[Join[StatesByRegionsEntity[[1]],StatesByRegionsEntity[[2]]]]},
		{GeoStyling["OutlineMap",Lighter[Blue]],Polygon[StatesByRegionsEntity[[3]]]},
		{GeoStyling["OutlineMap",Lighter[Blue,2/3]],Polygon[StatesByRegionsEntity[[4]]]},
		{GeoStyling["OutlineMap",Orange],Polygon[StatesByRegionsEntity[[5]]]},
		{GeoStyling["OutlineMap",Lighter[Orange]],Polygon[StatesByRegionsEntity[[6]]]},
		{GeoStyling["OutlineMap",Lighter[Orange,2/3]],Polygon[StatesByRegionsEntity[[7]]]},
		{GeoStyling["OutlineMap",Lighter[Red]],Polygon[Join[StatesByRegionsEntity[[8]],StatesByRegionsEntity[[9]]]]}
	};
	
	For[r=2,r<=8,r++,
		For[i=2,i<=8,i++,
			If[i!=r,
				AppendTo[styles,{Opacity[normalized[[r,i]]],Black,CapForm["Round"],Thickness[maxThickness normalized[[r,i]]],GeoPath[{regionHubEntity[[r]],regionHubEntity[[i]]}]}]
			]
		]
	];
	styles=Flatten[styles];
	
	Print[GeoGraphics[styles,GeoBackground->None]];
	Export["regionIBD_AfrAm-AfrAm_bin"<>ToString[binNumber]<>"_map.pdf",GeoGraphics[styles,GeoBackground->None]];
	
	(* region-specific relatedness *)
	For[r=2,r<=8,r++,
		minIBD=Min[Select[interRegionIBD[[r]],#!=0&]];
		maxIBD=Max[Select[interRegionIBD[[r]],#!=0&]];
		normalized=Map[If[#!=0,(1-cutoff) (#-minIBD)/(maxIBD-minIBD)+cutoff,0]&,interRegionIBD[[r]]];
		
		ClearAll[styles];
		styles=
		{
			{GeoStyling["OutlineMap",Darker[Blue]],Polygon[Join[StatesByRegionsEntity[[1]],StatesByRegionsEntity[[2]]]]},
			{GeoStyling["OutlineMap",Lighter[Blue]],Polygon[StatesByRegionsEntity[[3]]]},
			{GeoStyling["OutlineMap",Lighter[Blue,2/3]],Polygon[StatesByRegionsEntity[[4]]]},
			{GeoStyling["OutlineMap",Orange],Polygon[StatesByRegionsEntity[[5]]]},
			{GeoStyling["OutlineMap",Lighter[Orange]],Polygon[StatesByRegionsEntity[[6]]]},
			{GeoStyling["OutlineMap",Lighter[Orange,2/3]],Polygon[StatesByRegionsEntity[[7]]]},
			{GeoStyling["OutlineMap",Lighter[Red]],Polygon[Join[StatesByRegionsEntity[[8]],StatesByRegionsEntity[[9]]]]}
		};
		
		For[i=2,i<=8,i++,
			If[i!=r,
				AppendTo[styles,{Opacity[normalized[[i]]],Black,Thickness[maxThickness normalized[[i]]],CapForm["Round"],GeoPath[{regionHubEntity[[r]],regionHubEntity[[i]]}]}]
			]
		];
		styles=Flatten[styles];
		
		g[r]=GeoGraphics[styles,GeoBackground->None];
		Print[g[r]]
	]
,{binNumber,num}];

(*
For[r=2,r\[LessEqual]8,r++,
	Export["regionIBD_AfrAm-AfrAm_bin"<>ToString[binNumber]<>"_region"<>ToString[r]<>"_map.pdf",g[r]]
]
Remove[g]
*)



(* MAPS -- FOR EUROPEAN-AMERICANS -- RELATEDNESS BETWEEN REGIONS FOR DIFFERENT BINS -- USING IBD LENGTH *)
Print[Style["EurAm vs EurAm",18,Blue]]

cutoff=0.25; (* rescaling total IBD length to [cutoff, 1] *)
maxThickness=0.0125;

StatesByRegions=
{
	{"Connecticut","Maine","Massachusetts","New Hampshire","Rhode Island","Vermont"},
	{"New Jersey","New York","Pennsylvania"},
	{"Illinois","Indiana","Michigan","Ohio","Wisconsin"},
	{"Iowa","Kansas","Minnesota","Missouri","Nebraska","North Dakota","South Dakota"},
	{"Delaware","District Of Columbia","Florida","Georgia","Maryland","North Carolina","South Carolina","Virginia","West Virginia"},
	{"Alabama","Kentucky","Mississippi","Tennessee"},
	{"Arkansas","Louisiana","Oklahoma","Texas"},
	{"Arizona","Colorado","Idaho","Montana","Nevada","New Mexico","Utah","Wyoming"},
	{(*"Alaska",*)"California",(*"Hawaii",*)"Oregon","Washington"}
};
StatesByRegionsEntity=Table[Entity["AdministrativeDivision",{StringReplace[#," "->""],"UnitedStates"}]&/@StatesByRegions[[i]],{i,9}];

regionHub=
{
	{"Boston","Massachusetts","UnitedStates"},
	(*{"Philadelphia","Pennsylvania","UnitedStates"}*){"Burlington","Vermont","UnitedStates"},
	{"Chicago","Illinois","UnitedStates"},
	(*{"Omaha","Nebraska","UnitedStates"}*)(*{"Minneapolis","Minnesota","UnitedStates"}*){"Fargo","NorthDakota","UnitedStates"},
	(*{"Atlanta","Georgia","UnitedStates"}*){"Jacksonville","Florida","UnitedStates"},
	(*{"Jackson","Mississippi","UnitedStates"}*){"Waynesboro","Mississippi","UnitedStates"},
	(*{"NewOrleans","Louisiana","UnitedStates"}*){"Houston","Texas","UnitedStates"},
	(*{"Denver","Colorado","UnitedStates"}*){"SaltLakeCity","Utah","UnitedStates"},
	{"SanFrancisco","California","UnitedStates"}
};
regionHubEntity=Entity["City",#]&/@regionHub;

Table[
	Print[Style["\n\n\nBIN "<>ToString[binNumber]<>": "<>ToString[intervals[[binNumber]]]<>"cM",26,Red]];
	
	interRegionIBD=(EUEUIBDLength[[binNumber]]+Transpose[EUEUIBDLength[[binNumber]]]-2 DiagonalMatrix[Diagonal[EUEUIBDLength[[binNumber]]]]);
	minIBD=Min[Select[Flatten[interRegionIBD],#!=0&]];
	maxIBD=Max[interRegionIBD];
	(* fixing the scaling of the lines for AfrAm-AfrAm and EurAm-EurAm plots to the same scale *)
	(*interRegionIBDAfrAmAfrAm=(AfrAmAfrAmIBDLength\[LeftDoubleBracket]binNumber\[RightDoubleBracket]+Transpose[AfrAmAfrAmIBDLength\[LeftDoubleBracket]binNumber\[RightDoubleBracket]]-2 DiagonalMatrix[Diagonal[AfrAmAfrAmIBDLength\[LeftDoubleBracket]binNumber\[RightDoubleBracket]]]) maskAfrAmAfrAm(* apply the mask constructed above *);
	minIBD=Min[Min[Select[Flatten[interRegionIBD],#\[NotEqual]0&]],Min[Select[Flatten[interRegionIBDAfrAmAfrAm],#\[NotEqual]0&]]];
	maxIBD=Max[Max[interRegionIBD],Max[interRegionIBDAfrAmAfrAm]];*)
	
	normalized=Map[If[#!=0,(1-cutoff) (#-minIBD)/(maxIBD-minIBD)+cutoff,0]&,interRegionIBD,{2}];
	
	ClearAll[styles];
	styles=
	{
		{GeoStyling["OutlineMap",Darker[Blue]],Polygon[Join[StatesByRegionsEntity[[1]],StatesByRegionsEntity[[2]]]]},
		{GeoStyling["OutlineMap",Lighter[Blue]],Polygon[StatesByRegionsEntity[[3]]]},
		{GeoStyling["OutlineMap",Lighter[Blue,2/3]],Polygon[StatesByRegionsEntity[[4]]]},
		{GeoStyling["OutlineMap",Orange],Polygon[StatesByRegionsEntity[[5]]]},
		{GeoStyling["OutlineMap",Lighter[Orange]],Polygon[StatesByRegionsEntity[[6]]]},
		{GeoStyling["OutlineMap",Lighter[Orange,2/3]],Polygon[StatesByRegionsEntity[[7]]]},
		{GeoStyling["OutlineMap",Lighter[Red]],Polygon[Join[StatesByRegionsEntity[[8]],StatesByRegionsEntity[[9]]]]}
	};
	
	For[r=2,r<=8,r++,
		For[i=2,i<=8,i++,
			If[i!=r,
				AppendTo[styles,{Opacity[normalized[[r,i]]],Black,CapForm["Round"],Thickness[maxThickness normalized[[r,i]]],GeoPath[{regionHubEntity[[r]],regionHubEntity[[i]]}]}]
			]
		]
	];
	styles=Flatten[styles];
	
	Print[GeoGraphics[styles,GeoBackground->None]];
	Export["regionIBD_EurAm-EurAm_bin"<>ToString[binNumber]<>"_map.pdf",GeoGraphics[styles,GeoBackground->None]];
	
	(* region-specific relatedness *)
	For[r=2,r<=8,r++,
		minIBD=Min[Select[interRegionIBD[[r]],#!=0&]];
		maxIBD=Max[Select[interRegionIBD[[r]],#!=0&]];
		normalized=Map[If[#!=0,(1-cutoff) (#-minIBD)/(maxIBD-minIBD)+cutoff,0]&,interRegionIBD[[r]]];
		
		ClearAll[styles];
		styles=
		{
			{GeoStyling["OutlineMap",Darker[Blue]],Polygon[Join[StatesByRegionsEntity[[1]],StatesByRegionsEntity[[2]]]]},
			{GeoStyling["OutlineMap",Lighter[Blue]],Polygon[StatesByRegionsEntity[[3]]]},
			{GeoStyling["OutlineMap",Lighter[Blue,2/3]],Polygon[StatesByRegionsEntity[[4]]]},
			{GeoStyling["OutlineMap",Orange],Polygon[StatesByRegionsEntity[[5]]]},
			{GeoStyling["OutlineMap",Lighter[Orange]],Polygon[StatesByRegionsEntity[[6]]]},
			{GeoStyling["OutlineMap",Lighter[Orange,2/3]],Polygon[StatesByRegionsEntity[[7]]]},
			{GeoStyling["OutlineMap",Lighter[Red]],Polygon[Join[StatesByRegionsEntity[[8]],StatesByRegionsEntity[[9]]]]}
		};
		
		For[i=2,i<=8,i++,
			If[i!=r,
				AppendTo[styles,{Opacity[normalized[[i]]],Black,Thickness[maxThickness normalized[[i]]],CapForm["Round"],GeoPath[{regionHubEntity[[r]],regionHubEntity[[i]]}]}]
			]
		];
		styles=Flatten[styles];
		
		g[r]=GeoGraphics[styles,GeoBackground->None];
		Print[g[r]]
	]
,{binNumber,num}];

(*
For[r=2,r\[LessEqual]8,r++,
	Export["regionIBD_EurAm-EurAm_bin"<>ToString[binNumber]<>"_region"<>ToString[r]<>"_map.pdf",g[r]]
]
Remove[g]
*)






(* CENSUS DATA *)
(* DETAILS OF IPUMS DATA FIELDS *)
IPUMSdetails=DeleteCases[StringSplit[StringTrim[#]&/@Drop[Import["/PATH/TO/IPUMS/DATA/usa_00004.sps","List"],6],{"     ","    ","   ","  ","\""}],{}]

dictREGIONstart=Flatten[Position[IPUMSdetails,"/REGION"]][[1]]+1;
dictREGIONend=Flatten[Position[IPUMSdetails,"/STATEICP"]][[1]]-1;
dictREGIONlist=(#[[1]]->#[[-1]])&/@IPUMSdetails[[dictREGIONstart;;dictREGIONend]];
dictREGION=Association[dictREGIONlist]

dictSTATEstart=Flatten[Position[IPUMSdetails,"/STATEFIP"]][[1]]+1;
dictSTATEend=Flatten[Position[IPUMSdetails,"/CITY"]][[1]]-1;
dictSTATElist=(#[[1]]->#[[-1]])&/@IPUMSdetails[[dictSTATEstart;;dictSTATEend]];
dictSTATE=Association[dictSTATElist]

dictRACEstart=Flatten[Position[IPUMSdetails,"/RACE"]][[1]]+1;
dictRACEend=Flatten[Position[IPUMSdetails,"/RACED"]][[1]]-1;
dictRACElist=(#[[1]]->#[[-1]])&/@IPUMSdetails[[dictRACEstart;;dictRACEend]];
dictRACE=Association[dictRACElist]

dictBIRTHPLACEstart=Flatten[Position[IPUMSdetails,"/BPL"]][[1]]+1;
dictBIRTHPLACEend=Flatten[Position[IPUMSdetails,"/BPLD"]][[1]]-1;
dictBIRTHPLACElist=(#[[1]]->#[[-1]])&/@IPUMSdetails[[dictBIRTHPLACEstart;;dictBIRTHPLACEend]];
dictBIRTHPLACE=Association[dictBIRTHPLACElist]



(* LOCATION OF THE NECESSARY INFORMATION IN THE DATA STRING *)
YEARstart=1;
YEARend=4;
REGIONstart=25;
REGIONend=26;
STATEstart=29;
STATEend=30;
PERSONWEIGHTstart=40;
PERSONWEIGHTend=49;
AGEstart=50;
AGEend=52;
RACEstart=53;
RACEend=53;
BIRTHPLACEstart=57;
BIRTHPLACEend=59;



(* READING RAW CENSUS DATA *)
(* for large files, this gets super slow when it starts to use swap -- use the method below *)
(*censusData=ToString[#]&/@Import["/PATH/TO/IPUMS/DATA/usa_00003.dat.gz","List"];*)
Remove[stream,censusData,dataChunk]
stream=OpenRead["/PATH/TO/IPUMS/DATA/usa_00004.dat"];
censusData={};
While[dataChunk=!={},
	dataChunk=ReadList[stream,String,1000];
	AppendTo[censusData,dataChunk]
]
Close[stream];
Remove[stream,dataChunk]

censusData=Flatten[censusData];



(* EXTRACTING/REFORMATTING DETAILS FROM RAW CENSUS DATA *)
(* {YEAR, REGION, STATE, WEIGHT, AGE, RACE, BIRTHPLACE} *)
indivsData=
{
	ToExpression[StringTake[#,{YEARstart,YEARend}]],
	dictREGION[StringTake[#,{REGIONstart,REGIONend}]],
	dictSTATE[StringTake[#,{STATEstart,STATEend}]],
	ToExpression[StringTake[#,{PERSONWEIGHTstart,PERSONWEIGHTend}]]/100.0,
	ToExpression[StringTake[#,{AGEstart,AGEend}]],
	dictRACE[StringTake[#,{RACEstart,RACEend}]],
	dictBIRTHPLACE[StringTake[#,{BIRTHPLACEstart,BIRTHPLACEend}]]
}&/@censusData;



(* DICTIONARY FOR REGION NAME/NUMBER ASSOCIATION *)
dictREGIONname2number=Association[{"New England Division"->1,"Middle Atlantic Division"->2,"East North Central Div."->3,"West North Central Div."->4,"South Atlantic Division"->5,"East South Central Div."->6,"West South Central Div."->7,"Mountain Division"->8,"Pacific Division"->9}]
StatesByRegions=
{
	{"Connecticut","Maine","Massachusetts","New Hampshire","Rhode Island","Vermont"},
	{"New Jersey","New York","Pennsylvania"},
	{"Illinois","Indiana","Michigan","Ohio","Wisconsin"},
	{"Iowa","Kansas","Minnesota","Missouri","Nebraska","North Dakota","South Dakota"},
	{"Delaware","District of Columbia","Florida","Georgia","Maryland","North Carolina","South Carolina","Virginia","West Virginia"},
	{"Alabama","Kentucky","Mississippi","Tennessee"},
	{"Arkansas","Louisiana","Oklahoma","Texas"},
	{"Arizona","Colorado","Idaho","Montana","Nevada","New Mexico","Utah","Wyoming"},
	{(*"Alaska",*)"California",(*"Hawaii",*)"Oregon","Washington"}
};
dictSTATEtoREGION=Association[Flatten[Thread[#[[1]]->#[[2,1]]]&/@MapIndexed[List,StatesByRegions]]]



(* IDENTIFYING AFRICAN-AMERICANS AND EUROPEAN-AMERICANS IN THE CENSUS DATA *)
AfrAms=Select[indivsData,(#[[6]]=="Black/Negro")&]; (* labels defined by US Census Bureau *)
EurAms=Select[indivsData,(#[[6]]=="White")&]; (* labels defined by US Census Bureau *)



(* IDENTIFYING AFRICAN-AMERICANS AND EUROPEAN-AMERICANS WITHIN AGE RANGE 20-30 IN THE CENSUS DATA AND SORTING BY CENSUS DECADE *)
AfrAms\[LetterSpace]20\[LetterSpace]30=Select[AfrAms,(20<=#[[5]]<30)&];
AfrAms\[LetterSpace]20\[LetterSpace]30\[LetterSpace]censusDecades=Table[Select[AfrAms\[LetterSpace]20\[LetterSpace]30,(#[[1]]==1900+10 (i-1))&],{i,1,9}];

Length[AfrAms]
Length[AfrAms\[LetterSpace]20\[LetterSpace]30]
Length[#]&/@AfrAms\[LetterSpace]20\[LetterSpace]30\[LetterSpace]censusDecades

EurAms\[LetterSpace]20\[LetterSpace]30=Select[EurAms,(20<=#[[5]]<30)&];
EurAms\[LetterSpace]20\[LetterSpace]30\[LetterSpace]censusDecades=Table[Select[EurAms\[LetterSpace]20\[LetterSpace]30,(#[[1]]==1900+10 (i-1))&],{i,1,9}];

Length[EurAms]
Length[EurAms\[LetterSpace]20\[LetterSpace]30]
Length[#]&/@EurAms\[LetterSpace]20\[LetterSpace]30\[LetterSpace]censusDecades



(* IDENTIFYING MIGRATIONS FROM OUTSIDE OF THE CONTINENTAL US INTO THE US *)
notFromUS\[LetterSpace]indivsData=
Select[indivsData,!MemberQ[{"Alabama",(*"Alaska",*)"Arizona","Arkansas","California","Colorado","Connecticut","Delaware","District of Columbia","Florida","Georgia",(*"Hawaii",*)"Idaho","Illinois","Indiana","Iowa","Kansas","Kentucky","Louisiana","Maine","Maryland","Massachusetts","Michigan","Minnesota","Mississippi","Missouri","Montana","Nebraska","Nevada","New Hampshire","New Jersey","New Mexico","New York","North Carolina","North Dakota","Ohio","Oklahoma","Oregon","Pennsylvania","Rhode Island","South Carolina","South Dakota","Tennessee","Texas","Utah","Vermont","Virginia","Washington","West Virginia","Wisconsin","Wyoming"},#[[7]]]&];
Length[notFromUS\[LetterSpace]indivsData]

notFromUS\[LetterSpace]indivsData\[LetterSpace]AfrAm=Select[notFromUS\[LetterSpace]indivsData,#[[6]]=="Black/Negro"&]; (* labels defined by US Census Bureau *)
notFromUS\[LetterSpace]indivsData\[LetterSpace]EurAm=Select[notFromUS\[LetterSpace]indivsData,#[[6]]=="White"&]; (* labels defined by US Census Bureau *)
Length[notFromUS\[LetterSpace]indivsData\[LetterSpace]AfrAm]
Length[notFromUS\[LetterSpace]indivsData\[LetterSpace]EurAm]

Sort[Tally[notFromUS\[LetterSpace]indivsData\[LetterSpace]AfrAm[[All,3]]]]
Sort[Tally[notFromUS\[LetterSpace]indivsData\[LetterSpace]EurAm[[All,3]]]]
(* migrations to census regions, for AfrAms and EurAms *)
{dictSTATEtoREGION[#[[1]]],#[[2]]}&/@Tally[notFromUS\[LetterSpace]indivsData\[LetterSpace]AfrAm[[All,3]]]
Table[Total[Select[{dictSTATEtoREGION[#[[1]]],#[[2]]}&/@Tally[notFromUS\[LetterSpace]indivsData\[LetterSpace]AfrAm[[All,3]]],#[[1]]==regionNumber&][[All,2]]],{regionNumber,1,9}]
{dictSTATEtoREGION[#[[1]]],#[[2]]}&/@Tally[notFromUS\[LetterSpace]indivsData\[LetterSpace]EurAm[[All,3]]]
Table[Total[Select[{dictSTATEtoREGION[#[[1]]],#[[2]]}&/@Tally[notFromUS\[LetterSpace]indivsData\[LetterSpace]EurAm[[All,3]]],#[[1]]==regionNumber&][[All,2]]],{regionNumber,1,9}]



(* PROCESSING OUT_OF_US --> US MIGRATIONS, BASED ON AGE AND CENSUS YEAR *)
notFromUS\[LetterSpace]AfrAm\[LetterSpace]20\[LetterSpace]30=Select[notFromUS\[LetterSpace]indivsData\[LetterSpace]AfrAm,(20<=#[[5]]<30)&];
notFromUS\[LetterSpace]AfrAm\[LetterSpace]20\[LetterSpace]30\[LetterSpace]censusDecades=Table[Select[notFromUS\[LetterSpace]AfrAm\[LetterSpace]20\[LetterSpace]30,(#[[1]]==1900+10 (i-1))&],{i,1,9}];

notFromUS\[LetterSpace]EurAm\[LetterSpace]20\[LetterSpace]30=Select[notFromUS\[LetterSpace]indivsData\[LetterSpace]EurAm,(20<=#[[5]]<30)&];
notFromUS\[LetterSpace]EurAm\[LetterSpace]20\[LetterSpace]30\[LetterSpace]censusDecades=Table[Select[notFromUS\[LetterSpace]EurAm\[LetterSpace]20\[LetterSpace]30,(#[[1]]==1900+10 (i-1))&],{i,1,9}];

Length[notFromUS\[LetterSpace]AfrAm\[LetterSpace]20\[LetterSpace]30]
Tally[(dictREGIONname2number[#]&/@notFromUS\[LetterSpace]AfrAm\[LetterSpace]20\[LetterSpace]30[[All,2]])]
Tally[(dictSTATEtoREGION[#]&/@notFromUS\[LetterSpace]AfrAm\[LetterSpace]20\[LetterSpace]30[[All,3]])]
Print["\n"]
Length[notFromUS\[LetterSpace]EurAm\[LetterSpace]20\[LetterSpace]30]
Tally[(dictREGIONname2number[#]&/@notFromUS\[LetterSpace]EurAm\[LetterSpace]20\[LetterSpace]30[[All,2]])]
Tally[(dictSTATEtoREGION[#]&/@notFromUS\[LetterSpace]EurAm\[LetterSpace]20\[LetterSpace]30[[All,3]])]



(* IDENTIFYING MIGRATION PATTERNS IN CENSUS REGIONS *)
(* NOTE: using dictREGIONname2number[] to convert names of census regions to region numbers actually counts the people from Alaska and Hawaii as well, which is not what I do when using the IBD data -- DO NOT USE! *)
migrationsAfrAms=Table[0,{9},{9}];
migrationsEurAms=Table[0,{9},{9}];
For[i=1,i<=9,i++,
	For[j=1,j<=9,j++,
		Print[i," ",j];
		(* migrating from region 'i' to region 'j' *)
		(*migrationsAfrAms\[LeftDoubleBracket]i,j\[RightDoubleBracket]=Length[Select[#,(dictREGIONname2number[#\[LeftDoubleBracket]2\[RightDoubleBracket]]\[Equal]j&&dictSTATEtoREGION[#\[LeftDoubleBracket]7\[RightDoubleBracket]]\[Equal]i)&]]&/@AfrAms\[LetterSpace]20\[LetterSpace]30\[LetterSpace]censusDecades;*)
		(*migrationsAfrAms\[LeftDoubleBracket]i,j\[RightDoubleBracket]=Total[Select[#,(dictREGIONname2number[#\[LeftDoubleBracket]2\[RightDoubleBracket]]\[Equal]j&&dictSTATEtoREGION[#\[LeftDoubleBracket]7\[RightDoubleBracket]]\[Equal]i)&]\[LeftDoubleBracket]All,4\[RightDoubleBracket]]&/@AfrAms\[LetterSpace]20\[LetterSpace]30\[LetterSpace]censusDecades;*)(* see NOTE above *)
		migrationsAfrAms[[i,j]]=Total[Select[#,(dictSTATEtoREGION[#[[3]]]==j&&dictSTATEtoREGION[#[[7]]]==i)&][[All,4]]]&/@AfrAms\[LetterSpace]20\[LetterSpace]30\[LetterSpace]censusDecades;
		
		(*migrationsEurAms\[LeftDoubleBracket]i,j\[RightDoubleBracket]=Length[Select[#,(dictREGIONname2number[#\[LeftDoubleBracket]2\[RightDoubleBracket]]\[Equal]j&&dictSTATEtoREGION[#\[LeftDoubleBracket]7\[RightDoubleBracket]]\[Equal]i)&]]&/@EurAms\[LetterSpace]20\[LetterSpace]30\[LetterSpace]censusDecades*)
		(*migrationsEurAms\[LeftDoubleBracket]i,j\[RightDoubleBracket]=Total[Select[#,(dictREGIONname2number[#\[LeftDoubleBracket]2\[RightDoubleBracket]]\[Equal]j&&dictSTATEtoREGION[#\[LeftDoubleBracket]7\[RightDoubleBracket]]\[Equal]i)&]\[LeftDoubleBracket]All,4\[RightDoubleBracket]]&/@EurAms\[LetterSpace]20\[LetterSpace]30\[LetterSpace]censusDecades;*)(* see NOTE above *)
		migrationsEurAms[[i,j]]=Total[Select[#,(dictSTATEtoREGION[#[[3]]]==j&&dictSTATEtoREGION[#[[7]]]==i)&][[All,4]]]&/@EurAms\[LetterSpace]20\[LetterSpace]30\[LetterSpace]censusDecades
	]
]

notFromUS\[LetterSpace]migrationsAfrAms=Table[0,{9}];
notFromUS\[LetterSpace]migrationsEurAms=Table[0,{9}];
For[j=1,j<=9,j++,
	Print[j];
	(* migrating from outside of US to region 'j' *)
	notFromUS\[LetterSpace]migrationsAfrAms[[j]]=Total[Select[#,dictSTATEtoREGION[#[[3]]]==j&][[All,4]]]&/@notFromUS\[LetterSpace]AfrAm\[LetterSpace]20\[LetterSpace]30\[LetterSpace]censusDecades;
	notFromUS\[LetterSpace]migrationsEurAms[[j]]=Total[Select[#,dictSTATEtoREGION[#[[3]]]==j&][[All,4]]]&/@notFromUS\[LetterSpace]EurAm\[LetterSpace]20\[LetterSpace]30\[LetterSpace]censusDecades
]

Total[Total[notFromUS\[LetterSpace]migrationsAfrAms,{2}]/100(* samples have weight! *)]
Total[Total[notFromUS\[LetterSpace]migrationsEurAms,{2}]/100(* samples have weight! *)]
notFromUS\[LetterSpace]migrationsAfrAms//MatrixForm



(* CREATING THE COLLAPSED CENSUS REGIONS (1+2 --> 2 AND 8+9 --> 9) *)
migrationsCollapsedAfrAms=migrationsAfrAms;

(* migrations FROM regions 1 and 9... *)
migrationsCollapsedAfrAms[[2,2]]+=migrationsAfrAms[[1,1]];
migrationsCollapsedAfrAms[[2,{2,3,4,5,6,7,8}]]+=migrationsAfrAms[[1,{2,3,4,5,6,7,8}]];
migrationsCollapsedAfrAms[[2,8]]+=migrationsAfrAms[[1,9]];

migrationsCollapsedAfrAms[[8,2]]+=migrationsAfrAms[[9,1]];
migrationsCollapsedAfrAms[[8,{2,3,4,5,6,7,8}]]+=migrationsAfrAms[[9,{2,3,4,5,6,7,8}]];
migrationsCollapsedAfrAms[[8,8]]+=migrationsAfrAms[[9,9]];

migrationsCollapsedAfrAms[[1,All]]=Table[0,{9},{9}];
migrationsCollapsedAfrAms[[9,All]]=Table[0,{9},{9}];

(* migrations TO regions 1 and 9... *)
migrationsCollapsedAfrAms[[{2,3,4,5,6,7,8},2]]+=migrationsAfrAms[[{2,3,4,5,6,7,8},1]]; (* 1\[Rule]1 and 9\[Rule]1 migration is already taken care of above *)

migrationsCollapsedAfrAms[[{2,3,4,5,6,7,8},8]]+=migrationsAfrAms[[{2,3,4,5,6,7,8},9]]; (* 1\[Rule]9 and 9\[Rule]9 migration is already taken care of above *)

migrationsCollapsedAfrAms[[All,1]]=Table[0,{9},{9}];
migrationsCollapsedAfrAms[[All,9]]=Table[0,{9},{9}];

notFromUS\[LetterSpace]migrationsCollapsedAfrAms=notFromUS\[LetterSpace]migrationsAfrAms;
notFromUS\[LetterSpace]migrationsCollapsedAfrAms[[2]]+=notFromUS\[LetterSpace]migrationsAfrAms[[1]];
notFromUS\[LetterSpace]migrationsCollapsedAfrAms[[8]]+=notFromUS\[LetterSpace]migrationsAfrAms[[9]];
notFromUS\[LetterSpace]migrationsCollapsedAfrAms[[1]]=Table[0,{9}];
notFromUS\[LetterSpace]migrationsCollapsedAfrAms[[9]]=Table[0,{9}];

migrationsCollapsedEurAms=migrationsEurAms;

(* migrations FROM regions 1 and 9... *)
migrationsCollapsedEurAms[[2,2]]+=migrationsEurAms[[1,1]];
migrationsCollapsedEurAms[[2,{2,3,4,5,6,7,8}]]+=migrationsEurAms[[1,{2,3,4,5,6,7,8}]];
migrationsCollapsedEurAms[[2,8]]+=migrationsEurAms[[1,9]];

migrationsCollapsedEurAms[[8,2]]+=migrationsEurAms[[9,1]];
migrationsCollapsedEurAms[[8,{2,3,4,5,6,7,8}]]+=migrationsEurAms[[9,{2,3,4,5,6,7,8}]];
migrationsCollapsedEurAms[[8,8]]+=migrationsEurAms[[9,9]];

migrationsCollapsedEurAms[[1,All]]=Table[0,{9},{9}];
migrationsCollapsedEurAms[[9,All]]=Table[0,{9},{9}];

(* migrations TO regions 1 and 9... *)
migrationsCollapsedEurAms[[{2,3,4,5,6,7,8},2]]+=migrationsEurAms[[{2,3,4,5,6,7,8},1]]; (* 1\[Rule]1 and 9\[Rule]1 migration is already taken care of above *)

migrationsCollapsedEurAms[[{2,3,4,5,6,7,8},8]]+=migrationsEurAms[[{2,3,4,5,6,7,8},9]]; (* 1\[Rule]9 and 9\[Rule]9 migration is already taken care of above *)

migrationsCollapsedEurAms[[All,1]]=Table[0,{9},{9}];
migrationsCollapsedEurAms[[All,9]]=Table[0,{9},{9}];

notFromUS\[LetterSpace]migrationsCollapsedEurAms=notFromUS\[LetterSpace]migrationsEurAms;
notFromUS\[LetterSpace]migrationsCollapsedEurAms[[2]]+=notFromUS\[LetterSpace]migrationsEurAms[[1]];
notFromUS\[LetterSpace]migrationsCollapsedEurAms[[8]]+=notFromUS\[LetterSpace]migrationsEurAms[[9]];
notFromUS\[LetterSpace]migrationsCollapsedEurAms[[1]]=Table[0,{9}];
notFromUS\[LetterSpace]migrationsCollapsedEurAms[[9]]=Table[0,{9}];

(* an example *)
migrationsAfrAms[[5,2]]
ListPlot[Transpose[{Table[1900+10 (i-1),{i,1,9}],migrationsAfrAms[[5,2]]}],PlotRange->All]
Print["number of AfrAm migrants from 5 to 2 (in age range 20-30): ",Total[migrationsAfrAms[[5,2]]],"\n\n\n"]

migrationsAfrAms[[2,5]]
ListPlot[Transpose[{Table[1900+10 (i-1),{i,1,9}],migrationsAfrAms[[2,5]]}],PlotRange->All]
Print["number of AfrAm migrants from 2 to 5 (in age range 20-30): ",Total[migrationsAfrAms[[2,5]]],"\n\n\n"]

Print["ratio of the migrants from 2->5 to those from 5->2: ",N[Total[migrationsAfrAms[[2,5]]]/Total[migrationsAfrAms[[5,2]]]]]

Print["number of AfrAms who stayed in region 2 (in age range 20-30): ",migrationsAfrAms[[2,2]]]
Print["number of AfrAms who stayed in region 5 (in age range 20-30): ",migrationsAfrAms[[5,5]]]



(*
CALCULATING THE RELATEDNESS BETWEEN REGIONS
(in terms of proportion of genome of an individual in region i coming from an ancestor some generations back in region j -- see comments in the code below)
NOTE: the resulting matrix is not symmetric; migrations from 2->5 are different from 5->2, and this results in relatedness defined below to be directional.
NOTE: generations are defined as follows:
counter: 1 ==> start: 7, end: 9 (i.e., youngest, from 1960,1970,1980)
counter: 2 ==> start: 4, end: 6 (i.e., middle, from 1930,1940,1950)
counter: 3 ==> start: 1, end: 3 (i.e., oldest, from 1900,1910,1920)
*)
generationTime=30;
num\[LetterSpace]decades=Length[AfrAms\[LetterSpace]20\[LetterSpace]30\[LetterSpace]censusDecades];

counter=0;
For[end=num\[LetterSpace]decades,end>=1,end-=generationTime/10,counter++]

pAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest=Table[0,{counter},{9},{9}];
pAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]orig=Table[0,{counter},{9},{9}];
pCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest=Table[0,{counter},{9},{9}];
pCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]orig=Table[0,{counter},{9},{9}];

pEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest=Table[0,{counter},{9},{9}];
pEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]orig=Table[0,{counter},{9},{9}];
pCollapsedEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest=Table[0,{counter},{9},{9}];
pCollapsedEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]orig=Table[0,{counter},{9},{9}];

For[region\[LetterSpace]i=1,region\[LetterSpace]i<=9,region\[LetterSpace]i++,
	For[region\[LetterSpace]j=1,region\[LetterSpace]j<=9,region\[LetterSpace]j++,
		counter=0;
		For[end=num\[LetterSpace]decades,end>=1,end-=generationTime/10,
			counter++;
			If[end-(generationTime/10)+1<=0,start=1,start=end-(generationTime/10)+1];
			(*Print[start," ",end];*)
			
			pAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest[[counter,region\[LetterSpace]i,region\[LetterSpace]j]]=Total[migrationsAfrAms[[region\[LetterSpace]i,region\[LetterSpace]j,start;;end]]]/(Total[migrationsAfrAms[[All,region\[LetterSpace]j,start;;end]],2]+Total[notFromUS\[LetterSpace]migrationsAfrAms[[region\[LetterSpace]j,start;;end]]](* added for out-of-US\[Rule]j migrations *));
			pAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]orig[[counter,region\[LetterSpace]i,region\[LetterSpace]j]]=Total[migrationsAfrAms[[region\[LetterSpace]i,region\[LetterSpace]j,start;;end]]]/Total[migrationsAfrAms[[region\[LetterSpace]i,All,start;;end]],2](* no equivalent for the corresponding metric for offsprings *);
			
			pEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest[[counter,region\[LetterSpace]i,region\[LetterSpace]j]]=Total[migrationsEurAms[[region\[LetterSpace]i,region\[LetterSpace]j,start;;end]]]/(Total[migrationsEurAms[[All,region\[LetterSpace]j,start;;end]],2]+Total[notFromUS\[LetterSpace]migrationsEurAms[[region\[LetterSpace]j,start;;end]]](* added for out-of-US\[Rule]j migrations *));
			pEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]orig[[counter,region\[LetterSpace]i,region\[LetterSpace]j]]=Total[migrationsEurAms[[region\[LetterSpace]i,region\[LetterSpace]j,start;;end]]]/Total[migrationsEurAms[[region\[LetterSpace]i,All,start;;end]],2];
			
			pCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest[[counter,region\[LetterSpace]i,region\[LetterSpace]j]]=Total[migrationsCollapsedAfrAms[[region\[LetterSpace]i,region\[LetterSpace]j,start;;end]]]/If[Total[migrationsCollapsedAfrAms[[All,region\[LetterSpace]j,start;;end]],2]!=0,Total[migrationsCollapsedAfrAms[[All,region\[LetterSpace]j,start;;end]],2]+Total[notFromUS\[LetterSpace]migrationsCollapsedAfrAms[[region\[LetterSpace]j,start;;end]]](* added for out-of-US\[Rule]j migrations *),1];
			pCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]orig[[counter,region\[LetterSpace]i,region\[LetterSpace]j]]=Total[migrationsCollapsedAfrAms[[region\[LetterSpace]i,region\[LetterSpace]j,start;;end]]]/If[Total[migrationsCollapsedAfrAms[[region\[LetterSpace]i,All,start;;end]],2]!=0,Total[migrationsCollapsedAfrAms[[region\[LetterSpace]i,All,start;;end]],2],1];
			
			pCollapsedEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest[[counter,region\[LetterSpace]i,region\[LetterSpace]j]]=Total[migrationsCollapsedEurAms[[region\[LetterSpace]i,region\[LetterSpace]j,start;;end]]]/If[Total[migrationsCollapsedEurAms[[All,region\[LetterSpace]j,start;;end]],2]!=0,Total[migrationsCollapsedEurAms[[All,region\[LetterSpace]j,start;;end]],2]+Total[notFromUS\[LetterSpace]migrationsCollapsedEurAms[[region\[LetterSpace]j,start;;end]]](* added for out-of-US\[Rule]j migrations *),1];
			pCollapsedEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]orig[[counter,region\[LetterSpace]i,region\[LetterSpace]j]]=Total[migrationsCollapsedEurAms[[region\[LetterSpace]i,region\[LetterSpace]j,start;;end]]]/If[Total[migrationsCollapsedEurAms[[region\[LetterSpace]i,All,start;;end]],2]!=0,Total[migrationsCollapsedEurAms[[region\[LetterSpace]i,All,start;;end]],2],1];
		]
	]
]

normalizedMigrationsAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest=Apply[Dot,Table[pAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest[[i]],{i,counter,1,-1}]];
normalizedMigrationsAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]orig=Apply[Dot,Table[pAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]orig[[i]],{i,counter,1,-1}]];
normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest=Apply[Dot,Table[pCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest[[i]],{i,counter,1,-1}]];
normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]orig=Apply[Dot,Table[pCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]orig[[i]],{i,counter,1,-1}]];

normalizedMigrationsEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest=Apply[Dot,Table[pEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest[[i]],{i,counter,1,-1}]];
normalizedMigrationsEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]orig=Apply[Dot,Table[pEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]orig[[i]],{i,counter,1,-1}]];
normalizedMigrationsCollapsedEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest=Apply[Dot,Table[pCollapsedEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest[[i]],{i,counter,1,-1}]];
normalizedMigrationsCollapsedEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]orig=Apply[Dot,Table[pCollapsedEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]orig[[i]],{i,counter,1,-1}]];

symmetric\[LetterSpace]normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest=Table[0,{9},{9}];
symmetric\[LetterSpace]normalizedMigrationsCollapsedEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest=Table[0,{9},{9}];
For[i=2,i<=8,i++,
	For[j=2,j<=8,j++,
		symmetric\[LetterSpace]normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest[[i,j]]=Sum[1/migrationsCollapsedAfrAms[[k,k,1]] normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest[[k,i]] normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest[[k,j]],{k,2,8}];
		symmetric\[LetterSpace]normalizedMigrationsCollapsedEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest[[i,j]]=Sum[1/migrationsCollapsedEurAms[[k,k,1]] normalizedMigrationsCollapsedEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest[[k,i]] normalizedMigrationsCollapsedEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest[[k,j]],{k,2,8}]
	]
]



(* DIRECTIONAL RELATEDNESS BETWEEN REGIONS FOR AFRICAN-AMERICANS *)
ArrayPlot[normalizedMigrationsAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest,ColorFunctionScaling->False,Mesh->True] (* in the paper, we use "_dest" matrices *)
ArrayPlot[normalizedMigrationsAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]orig,ColorFunctionScaling->False,Mesh->True]
ArrayPlot[normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest,ColorFunctionScaling->False,Mesh->True] (* in the paper, we use "_dest" matrices *)
ArrayPlot[normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]orig,ColorFunctionScaling->False,Mesh->True]
ArrayPlot[symmetric\[LetterSpace]normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest,ColorFunctionScaling->False,Mesh->True]
ArrayPlot[UpperTriangularize[symmetric\[LetterSpace]normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest],ColorFunctionScaling->False,Mesh->True]

census\[LetterSpace]plot\[LetterSpace]dest=ArrayPlot[normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest[[{8,4,3,2,7,6,5},{8,4,3,2,7,6,5}]],Frame->True,FrameStyle->White,FrameTicks->{{None,MapIndexed[{First[#2],#1}&,regionByColor]},{None,MapIndexed[{First[#2],#1}&,regionByColor]}},ColorFunctionScaling->False,PlotRange->{Min[Select[Flatten[normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest],#!=0&]],Max[normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest]},PlotLegends->Placed[BarLegend[Automatic,LabelStyle->{Black,FontSize->30}],Bottom],ImageSize->750]
Export["regionCensus_AfrAm-AfrAm_ancestors_ReadByColumn_scale.pdf",census\[LetterSpace]plot\[LetterSpace]dest]

census\[LetterSpace]plot\[LetterSpace]orig=ArrayPlot[normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]orig[[{8,4,3,2,7,6,5},{8,4,3,2,7,6,5}]],Frame->True,FrameStyle->White,FrameTicks->{{None,MapIndexed[{First[#2],#1}&,regionByColor]},{None,MapIndexed[{First[#2],#1}&,regionByColor]}},ColorFunctionScaling->False,PlotRange->{Min[Select[Flatten[normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]orig],#!=0&]],Max[normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]orig]},PlotLegends->Placed[BarLegend[Automatic,LabelStyle->{Black,FontSize->30}],Bottom],ImageSize->750]
Export["regionCensus_AfrAm-AfrAm_descendants_ReadByRow_scale.pdf",census\[LetterSpace]plot\[LetterSpace]orig]

temp=ArrayPlot[Table[0,{7},{7}],Epilog->{MapIndexed[If[#1!=0,Text[Style[NumberForm[#1,{6,4}],26],Reverse[#2-1/2]]]&,Reverse[normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest[[{8,4,3,2,7,6,5},{8,4,3,2,7,6,5}]]],{2}]},Mesh->True,Frame->True,FrameStyle->White,FrameTicks->{{None,MapIndexed[{First[#2],#1}&,regionByColor]},{None,MapIndexed[{First[#2],#1}&,regionByColor]}},ImageSize->750]
Export["regionCensus_AfrAm-AfrAm_ancestors_ReadByColumn_number.pdf",temp]
Remove[temp]

temp=ArrayPlot[Table[0,{7},{7}],Epilog->{MapIndexed[If[#1!=0,Text[Style[NumberForm[#1,{6,3}],26],Reverse[#2-1/2]]]&,Reverse[normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]orig[[{8,4,3,2,7,6,5},{8,4,3,2,7,6,5}]]],{2}]},Mesh->True,Frame->True,FrameStyle->White,FrameTicks->{{None,MapIndexed[{First[#2],#1}&,regionByColor]},{None,MapIndexed[{First[#2],#1}&,regionByColor]}},ImageSize->750]
Export["regionCensus_AfrAm-AfrAm_descendants_ReadByRow_number.pdf",temp]
Remove[temp]

census\[LetterSpace]plot\[LetterSpace]dest=ArrayPlot[UpperTriangularize[symmetric\[LetterSpace]normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest[[{8,4,3,2,7,6,5},{8,4,3,2,7,6,5}]]],Frame->True,FrameStyle->White,FrameTicks->{{None,MapIndexed[{First[#2],#1}&,regionByColor]},{None,MapIndexed[{First[#2],#1}&,regionByColor]}},ColorFunctionScaling->False,PlotRange->{Min[Select[Flatten[symmetric\[LetterSpace]normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest],#!=0&]],Max[symmetric\[LetterSpace]normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest]},PlotLegends->Placed[BarLegend[Automatic,LabelStyle->{Black,FontSize->20}],Bottom],ImageSize->750]
Export["regionCensus_AfrAm-AfrAm_ancestors_symmetric_scale.pdf",census\[LetterSpace]plot\[LetterSpace]dest]

temp=ArrayPlot[Table[0,{7},{7}],Epilog->{MapIndexed[If[#1!=0,Text[Style[NumberForm[#1,{6,1}],20],Reverse[#2-1/2]]]&,Reverse[UpperTriangularize[symmetric\[LetterSpace]normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest[[{8,4,3,2,7,6,5},{8,4,3,2,7,6,5}]]]],{2}]},Mesh->True,Frame->True,FrameStyle->White,FrameTicks->{{None,MapIndexed[{First[#2],#1}&,regionByColor]},{None,MapIndexed[{First[#2],#1}&,regionByColor]}},ImageSize->750]
Export["regionCensus_AfrAm-AfrAm_ancestors_symmetric_number.pdf",temp]
Remove[temp]



(* RELATEDNESS MAPS BASED ON CENSUS DATA FOR AFRICAN-AMERICANS *)
Print[Style["AfrAm vs AfrAm",18,Blue]]

cutoff=0.25; (* rescaling total IBD length to [cutoff, 1] *)
maxThickness=0.0125;

StatesByRegions=
{
	{"Connecticut","Maine","Massachusetts","New Hampshire","Rhode Island","Vermont"},
	{"New Jersey","New York","Pennsylvania"},
	{"Illinois","Indiana","Michigan","Ohio","Wisconsin"},
	{"Iowa","Kansas","Minnesota","Missouri","Nebraska","North Dakota","South Dakota"},
	{"Delaware","District Of Columbia","Florida","Georgia","Maryland","North Carolina","South Carolina","Virginia","West Virginia"},
	{"Alabama","Kentucky","Mississippi","Tennessee"},
	{"Arkansas","Louisiana","Oklahoma","Texas"},
	{"Arizona","Colorado","Idaho","Montana","Nevada","New Mexico","Utah","Wyoming"},
	{(*"Alaska",*)"California",(*"Hawaii",*)"Oregon","Washington"}
};
StatesByRegionsEntity=Table[Entity["AdministrativeDivision",{StringReplace[#," "->""],"UnitedStates"}]&/@StatesByRegions[[i]],{i,9}];

regionHub=
{
	{"Boston","Massachusetts","UnitedStates"},
	(*{"Philadelphia","Pennsylvania","UnitedStates"}*){"Burlington","Vermont","UnitedStates"},
	{"Chicago","Illinois","UnitedStates"},
	(*{"Omaha","Nebraska","UnitedStates"}*)(*{"Minneapolis","Minnesota","UnitedStates"}*){"Fargo","NorthDakota","UnitedStates"},
	(*{"Atlanta","Georgia","UnitedStates"}*){"Jacksonville","Florida","UnitedStates"},
	(*{"Jackson","Mississippi","UnitedStates"}*){"Waynesboro","Mississippi","UnitedStates"},
	(*{"NewOrleans","Louisiana","UnitedStates"}*){"Houston","Texas","UnitedStates"},
	(*{"Denver","Colorado","UnitedStates"}*){"SaltLakeCity","Utah","UnitedStates"},
	{"SanFrancisco","California","UnitedStates"}
};
regionHubEntity=Entity["City",#]&/@regionHub;

interRegionIBD=(symmetric\[LetterSpace]normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest-DiagonalMatrix[Diagonal[symmetric\[LetterSpace]normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest]]) maskAfrAmAfrAm; (* apply the same mask that was used for plotting IBD on map for African-Americans *)
minIBD=Min[Select[Flatten[interRegionIBD],#!=0&]];
maxIBD=Max[interRegionIBD];
Print["min: ",minIBD]
Print["max: ",maxIBD]
normalized=Map[If[#!=0,(1-cutoff) (#-minIBD)/(maxIBD-minIBD)+cutoff,0]&,interRegionIBD,{2}];

ClearAll[styles];
styles=
{
	{GeoStyling["OutlineMap",Darker[Blue]],Polygon[Join[StatesByRegionsEntity[[1]],StatesByRegionsEntity[[2]]]]},
	{GeoStyling["OutlineMap",Lighter[Blue]],Polygon[StatesByRegionsEntity[[3]]]},
	{GeoStyling["OutlineMap",Lighter[Blue,2/3]],Polygon[StatesByRegionsEntity[[4]]]},
	{GeoStyling["OutlineMap",Orange],Polygon[StatesByRegionsEntity[[5]]]},
	{GeoStyling["OutlineMap",Lighter[Orange]],Polygon[StatesByRegionsEntity[[6]]]},
	{GeoStyling["OutlineMap",Lighter[Orange,2/3]],Polygon[StatesByRegionsEntity[[7]]]},
	{GeoStyling["OutlineMap",Lighter[Red]],Polygon[Join[StatesByRegionsEntity[[8]],StatesByRegionsEntity[[9]]]]}
};

For[r=2,r<=8,r++,
	For[i=2,i<=8,i++,
		If[i!=r,
			AppendTo[styles,{Opacity[normalized[[r,i]]],Black,CapForm["Round"],Thickness[maxThickness normalized[[r,i]]],GeoPath[{regionHubEntity[[r]],regionHubEntity[[i]]}]}]
		]
	]
];
styles=Flatten[styles];

Print[GeoGraphics[styles,GeoBackground->None]];
Export["regionCensus_AfrAm-AfrAm_ancestors_symmetric_map.pdf",GeoGraphics[styles,GeoBackground->None]];

(* region-specific relatedness *)
For[r=2,r<=8,r++,
	minIBD=Min[Select[interRegionIBD[[r]],#!=0&]];
	maxIBD=Max[Select[interRegionIBD[[r]],#!=0&]];
	normalized=Map[If[#!=0,(1-cutoff) (#-minIBD)/(maxIBD-minIBD)+cutoff,0]&,interRegionIBD[[r]]];
	
	ClearAll[styles];
	styles=
	{
		{GeoStyling["OutlineMap",Darker[Blue]],Polygon[Join[StatesByRegionsEntity[[1]],StatesByRegionsEntity[[2]]]]},
		{GeoStyling["OutlineMap",Lighter[Blue]],Polygon[StatesByRegionsEntity[[3]]]},
		{GeoStyling["OutlineMap",Lighter[Blue,2/3]],Polygon[StatesByRegionsEntity[[4]]]},
		{GeoStyling["OutlineMap",Orange],Polygon[StatesByRegionsEntity[[5]]]},
		{GeoStyling["OutlineMap",Lighter[Orange]],Polygon[StatesByRegionsEntity[[6]]]},
		{GeoStyling["OutlineMap",Lighter[Orange,2/3]],Polygon[StatesByRegionsEntity[[7]]]},
		{GeoStyling["OutlineMap",Lighter[Red]],Polygon[Join[StatesByRegionsEntity[[8]],StatesByRegionsEntity[[9]]]]}
	};
	
	For[i=2,i<=8,i++,
		If[i!=r,
			AppendTo[styles,{Opacity[normalized[[i]]],Black,Thickness[maxThickness normalized[[i]]],CapForm["Round"],GeoPath[{regionHubEntity[[r]],regionHubEntity[[i]]}]}]
		]
	];
	styles=Flatten[styles];
	
	g[r]=GeoGraphics[styles,GeoBackground->None];
	Print[g[r]]
]

(*
For[r=2,r\[LessEqual]8,r++,
	Export["regionCensus_AfrAm-AfrAm_ancestors_symmetric_region"<>ToString[r]<>"_map.pdf",g[r]]
]
Remove[g]
*)



(* DIRECTIONAL RELATEDNESS BETWEEN REGIONS FOR EUROPEAN-AMERICANS *)
ArrayPlot[normalizedMigrationsEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest,ColorFunctionScaling->False,Mesh->True] (* in the paper, we use "_dest" matrices *)
ArrayPlot[normalizedMigrationsEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]orig,ColorFunctionScaling->False,Mesh->True]
ArrayPlot[normalizedMigrationsCollapsedEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest,ColorFunctionScaling->False,Mesh->True] (* in the paper, we use "_dest" matrices *)
ArrayPlot[normalizedMigrationsCollapsedEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]orig,ColorFunctionScaling->False,Mesh->True]
ArrayPlot[symmetric\[LetterSpace]normalizedMigrationsCollapsedEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest,ColorFunctionScaling->False,Mesh->True]
ArrayPlot[UpperTriangularize[symmetric\[LetterSpace]normalizedMigrationsCollapsedEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest],ColorFunctionScaling->False,Mesh->True]

census\[LetterSpace]plot\[LetterSpace]dest=ArrayPlot[normalizedMigrationsCollapsedEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest[[{8,4,3,2,7,6,5},{8,4,3,2,7,6,5}]],Frame->True,FrameStyle->White,FrameTicks->{{None,MapIndexed[{First[#2],#1}&,regionByColor]},{None,MapIndexed[{First[#2],#1}&,regionByColor]}},ColorFunctionScaling->False,PlotRange->{Min[Select[Flatten[normalizedMigrationsCollapsedEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest],#!=0&]],Max[normalizedMigrationsCollapsedEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest]},PlotLegends->Placed[BarLegend[Automatic,LabelStyle->{Black,FontSize->30}],Bottom],ImageSize->750]
Export["regionCensus_EurAm-EurAm_ancestors_ReadByColumn_scale.pdf",census\[LetterSpace]plot\[LetterSpace]dest]

census\[LetterSpace]plot\[LetterSpace]orig=ArrayPlot[normalizedMigrationsCollapsedEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]orig[[{8,4,3,2,7,6,5},{8,4,3,2,7,6,5}]],Frame->True,FrameStyle->White,FrameTicks->{{None,MapIndexed[{First[#2],#1}&,regionByColor]},{None,MapIndexed[{First[#2],#1}&,regionByColor]}},ColorFunctionScaling->False,PlotRange->{Min[Select[Flatten[normalizedMigrationsCollapsedEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]orig],#!=0&]],Max[normalizedMigrationsCollapsedEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]orig]},PlotLegends->Placed[BarLegend[Automatic,LabelStyle->{Black,FontSize->30}],Bottom],ImageSize->750]
Export["regionCensus_EurAm-EurAm_descendants_ReadByRow_scale.pdf",census\[LetterSpace]plot\[LetterSpace]orig]

temp=ArrayPlot[Table[0,{7},{7}],Epilog->{MapIndexed[If[#1!=0,Text[Style[NumberForm[#1,{6,4}],26],Reverse[#2-1/2]]]&,Reverse[normalizedMigrationsCollapsedEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest[[{8,4,3,2,7,6,5},{8,4,3,2,7,6,5}]]],{2}]},Mesh->True,Frame->True,FrameStyle->White,FrameTicks->{{None,MapIndexed[{First[#2],#1}&,regionByColor]},{None,MapIndexed[{First[#2],#1}&,regionByColor]}},ImageSize->750]
Export["regionCensus_EurAm-EurAm_ancestors_ReadByColumn_number.pdf",temp]
Remove[temp]

temp=ArrayPlot[Table[0,{7},{7}],Epilog->{MapIndexed[If[#1!=0,Text[Style[NumberForm[#1,{6,3}],26],Reverse[#2-1/2]]]&,Reverse[normalizedMigrationsCollapsedEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]orig[[{8,4,3,2,7,6,5},{8,4,3,2,7,6,5}]]],{2}]},Mesh->True,Frame->True,FrameStyle->White,FrameTicks->{{None,MapIndexed[{First[#2],#1}&,regionByColor]},{None,MapIndexed[{First[#2],#1}&,regionByColor]}},ImageSize->750]
Export["regionCensus_EurAm-EurAm_descendants_ReadByRow_number.pdf",temp]
Remove[temp]

census\[LetterSpace]plot\[LetterSpace]dest=ArrayPlot[UpperTriangularize[symmetric\[LetterSpace]normalizedMigrationsCollapsedEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest[[{8,4,3,2,7,6,5},{8,4,3,2,7,6,5}]]],Frame->True,FrameStyle->White,FrameTicks->{{None,MapIndexed[{First[#2],#1}&,regionByColor]},{None,MapIndexed[{First[#2],#1}&,regionByColor]}},ColorFunctionScaling->False,PlotRange->{Min[Select[Flatten[symmetric\[LetterSpace]normalizedMigrationsCollapsedEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest],#!=0&]],Max[symmetric\[LetterSpace]normalizedMigrationsCollapsedEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest]},PlotLegends->Placed[BarLegend[Automatic,LabelStyle->{Black,FontSize->20}],Bottom],ImageSize->750]
Export["regionCensus_EurAm-EurAm_ancestors_symmetric_scale.pdf",census\[LetterSpace]plot\[LetterSpace]dest]

temp=ArrayPlot[Table[0,{7},{7}],Epilog->{MapIndexed[If[#1!=0,Text[Style[NumberForm[#1,{6,2}],19],Reverse[#2-1/2]]]&,Reverse[UpperTriangularize[symmetric\[LetterSpace]normalizedMigrationsCollapsedEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest[[{8,4,3,2,7,6,5},{8,4,3,2,7,6,5}]]]],{2}]},Mesh->True,Frame->True,FrameStyle->White,FrameTicks->{{None,MapIndexed[{First[#2],#1}&,regionByColor]},{None,MapIndexed[{First[#2],#1}&,regionByColor]}},ImageSize->750]
Export["regionCensus_EurAm-EurAm_ancestors_symmetric_number.pdf",temp]
Remove[temp]



(* RELATEDNESS MAPS BASED ON CENSUS DATA FOR EUROPEAN-AMERICANS *)
Print[Style["EurAm vs EurAm",18,Blue]]

cutoff=0.25; (* rescaling total IBD length to [cutoff, 1] *)
maxThickness=0.0125;

StatesByRegions=
{
	{"Connecticut","Maine","Massachusetts","New Hampshire","Rhode Island","Vermont"},
	{"New Jersey","New York","Pennsylvania"},
	{"Illinois","Indiana","Michigan","Ohio","Wisconsin"},
	{"Iowa","Kansas","Minnesota","Missouri","Nebraska","North Dakota","South Dakota"},
	{"Delaware","District Of Columbia","Florida","Georgia","Maryland","North Carolina","South Carolina","Virginia","West Virginia"},
	{"Alabama","Kentucky","Mississippi","Tennessee"},
	{"Arkansas","Louisiana","Oklahoma","Texas"},
	{"Arizona","Colorado","Idaho","Montana","Nevada","New Mexico","Utah","Wyoming"},
	{(*"Alaska",*)"California",(*"Hawaii",*)"Oregon","Washington"}
};
StatesByRegionsEntity=Table[Entity["AdministrativeDivision",{StringReplace[#," "->""],"UnitedStates"}]&/@StatesByRegions[[i]],{i,9}];

regionHub=
{
	{"Boston","Massachusetts","UnitedStates"},
	(*{"Philadelphia","Pennsylvania","UnitedStates"}*){"Burlington","Vermont","UnitedStates"},
	{"Chicago","Illinois","UnitedStates"},
	(*{"Omaha","Nebraska","UnitedStates"}*)(*{"Minneapolis","Minnesota","UnitedStates"}*){"Fargo","NorthDakota","UnitedStates"},
	(*{"Atlanta","Georgia","UnitedStates"}*){"Jacksonville","Florida","UnitedStates"},
	(*{"Jackson","Mississippi","UnitedStates"}*){"Waynesboro","Mississippi","UnitedStates"},
	(*{"NewOrleans","Louisiana","UnitedStates"}*){"Houston","Texas","UnitedStates"},
	(*{"Denver","Colorado","UnitedStates"}*){"SaltLakeCity","Utah","UnitedStates"},
	{"SanFrancisco","California","UnitedStates"}
};
regionHubEntity=Entity["City",#]&/@regionHub;

interRegionIBD=symmetric\[LetterSpace]normalizedMigrationsCollapsedEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest-DiagonalMatrix[Diagonal[symmetric\[LetterSpace]normalizedMigrationsCollapsedEurAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest]];
minIBD=Min[Select[Flatten[interRegionIBD],#!=0&]];
maxIBD=Max[interRegionIBD];
Print["min: ",minIBD]
Print["max: ",maxIBD]
normalized=Map[If[#!=0,(1-cutoff) (#-minIBD)/(maxIBD-minIBD)+cutoff,0]&,interRegionIBD,{2}];

ClearAll[styles];
styles=
{
	{GeoStyling["OutlineMap",Darker[Blue]],Polygon[Join[StatesByRegionsEntity[[1]],StatesByRegionsEntity[[2]]]]},
	{GeoStyling["OutlineMap",Lighter[Blue]],Polygon[StatesByRegionsEntity[[3]]]},
	{GeoStyling["OutlineMap",Lighter[Blue,2/3]],Polygon[StatesByRegionsEntity[[4]]]},
	{GeoStyling["OutlineMap",Orange],Polygon[StatesByRegionsEntity[[5]]]},
	{GeoStyling["OutlineMap",Lighter[Orange]],Polygon[StatesByRegionsEntity[[6]]]},
	{GeoStyling["OutlineMap",Lighter[Orange,2/3]],Polygon[StatesByRegionsEntity[[7]]]},
	{GeoStyling["OutlineMap",Lighter[Red]],Polygon[Join[StatesByRegionsEntity[[8]],StatesByRegionsEntity[[9]]]]}
};

For[r=2,r<=8,r++,
	For[i=2,i<=8,i++,
		If[i!=r,
			AppendTo[styles,{Opacity[normalized[[r,i]]],Black,CapForm["Round"],Thickness[maxThickness normalized[[r,i]]],GeoPath[{regionHubEntity[[r]],regionHubEntity[[i]]}]}]
		]
	]
];
styles=Flatten[styles];

Print[GeoGraphics[styles,GeoBackground->None]];
Export["regionCensus_EurAm-EurAm_ancestors_symmetric_map.pdf",GeoGraphics[styles,GeoBackground->None]];

(* region-specific relatedness *)
For[r=2,r<=8,r++,
	minIBD=Min[Select[interRegionIBD[[r]],#!=0&]];
	maxIBD=Max[Select[interRegionIBD[[r]],#!=0&]];
	normalized=Map[If[#!=0,(1-cutoff) (#-minIBD)/(maxIBD-minIBD)+cutoff,0]&,interRegionIBD[[r]]];
	
	ClearAll[styles];
	styles=
	{
		{GeoStyling["OutlineMap",Darker[Blue]],Polygon[Join[StatesByRegionsEntity[[1]],StatesByRegionsEntity[[2]]]]},
		{GeoStyling["OutlineMap",Lighter[Blue]],Polygon[StatesByRegionsEntity[[3]]]},
		{GeoStyling["OutlineMap",Lighter[Blue,2/3]],Polygon[StatesByRegionsEntity[[4]]]},
		{GeoStyling["OutlineMap",Orange],Polygon[StatesByRegionsEntity[[5]]]},
		{GeoStyling["OutlineMap",Lighter[Orange]],Polygon[StatesByRegionsEntity[[6]]]},
		{GeoStyling["OutlineMap",Lighter[Orange,2/3]],Polygon[StatesByRegionsEntity[[7]]]},
		{GeoStyling["OutlineMap",Lighter[Red]],Polygon[Join[StatesByRegionsEntity[[8]],StatesByRegionsEntity[[9]]]]}
	};
	
	For[i=2,i<=8,i++,
		If[i!=r,
			AppendTo[styles,{Opacity[normalized[[i]]],Black,Thickness[maxThickness normalized[[i]]],CapForm["Round"],GeoPath[{regionHubEntity[[r]],regionHubEntity[[i]]}]}]
		]
	];
	styles=Flatten[styles];
	
	g[r]=GeoGraphics[styles,GeoBackground->None];
	Print[g[r]]
]

(*
For[r=2,r\[LessEqual]8,r++,
	Export["regionCensus_EurAm-EurAm_ancestors_symmetric_region"<>ToString[r]<>"_map.pdf",g[r]]
]
Remove[g]
*)






(* CORRELATION BETWEEN GENOMIC AND CENSUS DATA *)
(* MANTEL TEST *)
correlationPearson[m1_,m2_]:=
Module[{m1Vect,m2Vect,nElements},
	m1Vect=Flatten[m1];
	m2Vect=Flatten[m2];
	nElements=Length[m1Vect];
	
	Return[N[(nElements m1Vect.m2Vect-Total[m1Vect] Total[m2Vect])/Sqrt[(nElements m1Vect.m1Vect-Total[m1Vect]^2) (nElements m2Vect.m2Vect-Total[m2Vect]^2)]]]
]

MantelTest[m1_?MatrixQ,m2_?MatrixQ,nPerm_?IntegerQ,keepDiag_?BooleanQ]:=
Module[{i,scoreOrig,scorePerm,nSuccess,dim,dimVect,permutation,m2Vect},
	scoreOrig=correlationPearson[m1,m2];
	
	dim=Length[m1];
	
	nSuccess=0;
	For[i=1,i<=nPerm,i++,
		If[TrueQ[keepDiag],
			permutation=RandomSample[Range[dim],dim];
			scorePerm=correlationPearson[m1,m2[[permutation,permutation]]]
		,
			m2Vect=Flatten[m2];
			dimVect=Length[m2Vect];
			permutation=RandomSample[Range[dimVect],dimVect];
			scorePerm=correlationPearson[m1,Partition[m2Vect[[permutation]],dim,dim]]
		];
		
		If[scoreOrig>=scorePerm,nSuccess++]
	];
	
	Return[{scoreOrig,nSuccess/nPerm,N[1-nSuccess/nPerm]}]
]



(* CALCULATING CORRELATIONS FOR THE 3x3 SOUTH-to-NORTH IBD MATRIX (genome data versus census data) *)
m1=AfrAmAfrAmIBDLength[[3]][[2;;4,5;;7]];
m2=symmetric\[LetterSpace]normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest[[2;;4,5;;7]];
ArrayPlot[m1,ColorFunctionScaling->False]
ArrayPlot[m2,ColorFunctionScaling->False]
MantelTest[m1,m2,9!,False]



(* CALCULATING CORRELATIONS FOR THE 4x3 SOUTH-to-(NORTH--and--WEST) IBD MATRIX (genome data versus census data) *)
m1=Join[AfrAmAfrAmIBDLength[[3]][[2;;4,5;;7]],{AfrAmAfrAmIBDLength[[3]][[5;;7,8]]}];
m2=symmetric\[LetterSpace]normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest[[{2,3,4,8},5;;7]];
ArrayPlot[m1,ColorFunctionScaling->False]
ArrayPlot[m2,ColorFunctionScaling->False]

MantelTest[m1,m2,10^7(*12! not feasible*),False]



(* CALCULATING CORRELATIONS FOR THE 3x3 SOUTH-to-NORTH (no WEST included) IBD MATRIX (genomic data using birth region versus the third (i.e., oldest -- 1900,1910,1920) generation relatedness matrix from census) *)
(* data from the AfrAmAfrAmIBDLength\[LeftDoubleBracket]3\[RightDoubleBracket]\[LeftDoubleBracket]2;;4,5;;7\[RightDoubleBracket] in ..._IBDbins.backup_20150309_v2.nb file *)
temp1={{0.04771368589682744`,0.03717620806950635`,0.030179589907955294`},{0.03670579615014729`,0.04874003181622716`,0.03992937873303168`},{0.031891161037852785`,0.055071169764833935`,0.04876065335716635`}};
normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest\[LetterSpace]gen3=pCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest[[3]];
symmetric\[LetterSpace]normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest\[LetterSpace]gen3=Table[0,{9},{9}];
For[i=2,i<=8,i++,
	For[j=2,j<=8,j++,
		symmetric\[LetterSpace]normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest\[LetterSpace]gen3[[i,j]]=Sum[1/migrationsCollapsedAfrAms[[k,k,1]] normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest\[LetterSpace]gen3[[k,i]] normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest\[LetterSpace]gen3[[k,j]],{k,2,8}]
	]
]
temp2=symmetric\[LetterSpace]normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest\[LetterSpace]gen3[[2;;4,5;;7]];
ArrayPlot[temp1,ColorFunctionScaling->False]
ArrayPlot[temp2,ColorFunctionScaling->False]

MantelTest[temp1,temp2,9!,False]



(* CALCULATING CORRELATIONS FOR THE 4x3 SOUTH-to-(NORTH--and--WEST) IBD MATRIX (genomic data using birth region versus the third (i.e., oldest -- 1900,1910,1920) generation relatedness matrix from census) *)
(* data from the Join[AfrAmAfrAmIBDLength\[LeftDoubleBracket]3\[RightDoubleBracket]\[LeftDoubleBracket]2;;4,5;;7\[RightDoubleBracket],{AfrAmAfrAmIBDLength\[LeftDoubleBracket]3\[RightDoubleBracket]\[LeftDoubleBracket]5;;7,8\[RightDoubleBracket]}] in ..._IBDbins.backup_20150309_v2.nb file *)
temp1={{0.04771368589682744`,0.03717620806950635`,0.030179589907955294`},{0.03670579615014729`,0.04874003181622716`,0.03992937873303168`},{0.031891161037852785`,0.055071169764833935`,0.04876065335716635`},{0.02576722113254151`,0.03447599531027599`,0.04050492016530479`}};
normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest\[LetterSpace]gen3=pCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest[[3]];
symmetric\[LetterSpace]normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest\[LetterSpace]gen3=Table[0,{9},{9}];
For[i=2,i<=8,i++,
	For[j=2,j<=8,j++,
		symmetric\[LetterSpace]normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest\[LetterSpace]gen3[[i,j]]=Sum[1/migrationsCollapsedAfrAms[[k,k,1]] normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest\[LetterSpace]gen3[[k,i]] normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest\[LetterSpace]gen3[[k,j]],{k,2,8}]
	]
]
temp2=symmetric\[LetterSpace]normalizedMigrationsCollapsedAfrAms\[LetterSpace]MATRIXFORM\[LetterSpace]dest\[LetterSpace]gen3[[{2,3,4,8},5;;7]];
ArrayPlot[temp1,ColorFunctionScaling->False]
ArrayPlot[temp2,ColorFunctionScaling->False]

MantelTest[temp1,temp2,10^7(*12! not feasible*),False]