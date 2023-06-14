# vaRHC
R package to automate, as far as possible, the variant classification process for hereditary cancer genes. 


## Introduction
Variant classification is a manual complex long process that combines information of distinct nature. An accurate classification is necessary to ensure proper genetic counselling and personalised risk estimation. 
In 2015, the American College of Molecular Genetics and Genomics (ACMG) together with the Association of Molecular Pathologists published generic guidelines to standardize and provide an objective framework to evaluate variant pathogenicity in Mendelian disease.
Later, specific guidelines have been published for some genes by collaborative groups. 

`vaRHC` has been developed to automate as much as possible the process of variant classification in hereditary cancer.
It follows gene-specific guidelines for *APC*, *ATM*, *CDH1*, *CHEK2*, *MLH1*, *MSH2*, *MSH6*, *PALB2*, *PMS2*, *PTEN* and *TP53*, and the updated general ACMG/AMP rules for other cancer susceptibility genes. The final classification is obtained according to Tavtigian’s natural scoring Bayesian-based metastructure (Tavtigian et al., 2020) but also considering some of the proposed CanVIG-UK incompatibilities.

The current version of the package works for single substitutions, deletions and insertions up to 25 bp, intronic variants and 5’ or 3’-UTR variants 25 bp.

## How to use it
Installation and user guide are available at the [vignette](https://htmlpreview.github.io/?https://github.com/emunte/vaRHC/blob/main/inst/doc/vaRHC.html).

## Citation
vaRHC was developed by Elisabet Munté Roca. The paper is now being submitted. 

## Legal advice and privacy policy for users
### Legal advice
#### Terms and conditions 
Access and navigation through the web application ascribe the User condition and imply the full and unreserved acceptance of each of the provisions included in this Legal Warning and Privacy Policy.  
Also, for any doubt the user has regarding using the software, he can contact emunte@idibell.cat. 

#### User's rights and obligations 
The user has the right to browse the website, observe the rules established in its notices and policies at all times, and the right to have your information processed with respect and maximum privacy, which is why we have provided adequate Privacy and Data Protection Policy. <br>

The user must always respect the terms and conditions established in this legal notice. Expressly, the user states he will use the portal diligently and assume any liability that may arise from the breach of rules. <bR>

Users will be obliged to make reasonable use of services or content, under the principle of good faith and respect for the legality in force, morality, public order, good customs, third-party rights, or IDIBELL itself, all according to the possibilities and objectives for which they are conceived. <br>

IDIBELL does not assume direct or indirect responsibilities for emerging damage or enduring profit resulting from users or third parties' misuse of services or content. <br>

The user undertakes to indemnify and hold harmless the website for any damage, prejudice, penalty, fine, penalty, or compensation that may have to have the website. 

### Confidentiality rules
Your data may be collected when the user browses the software, and, in this case, IDIBELL will respect the provisions of the Privacy Policy. <br>

The application can automatically detect the user's IP address and domain name. An IP address is a number automatically assigned to a computer when it connects to the Internet. All this information is registered in a server activity file duly registered that allows the subsequent processing of the data to obtain only statistical measurements that allow knowing the number of page impressions, the number of visits made to the web servers, the order of visits, the access point, and to assess other Website' aspects. 

##### Software of the web and its design 
IDIBELL is the owner of this software because it has been developed by its hired staff. 

##### Responsibility

IDIBELL does not guarantee the continuous and permanent availability of the services, and IDIBELL has not any responsibility for possible damages caused as a result of the lack of availability of the service due to force majeure or errors in the telematic networks of transfer data, alien to their will, or disconnections made for improvement or maintenance of equipment and computer systems.<br>

Also, not responsible for possible omissions, loss of data, settings, improper access, or breach of confidentiality that have their origin in technical problems, communications, or human failures, caused by third parties or not attributable to the software. Nor will it be liable for damages caused by computer attacks or caused by viruses that affect computer programs, communications systems, or equipment used by the website but manufactured or provided by a third party. <br>

IDIBELL will not be able to be responsible for any service suspension that is determined by a internal decision. <br> 

IDIBELL will be able to its discretion, deny, withdraw, suspend and/or block at any time, and without warning, access to information and services to users who breach the present rules except where the Act expressly imposes the opposite, and exclusively with the measure and extent to which it assesses, no liability is guaranteed or assumed for damages caused by software, data and portal services. <br>

##### Applicable law and jurisdiction

IDIBELL also reserved the right to file civil or criminal actions it deems appropriate for the improper use of its website and content or breaking these conditions.<br>

Current Spanish regulations will govern the relationships between IDIBELL and User, and any dispute will be submitted to the Court and Tribunals of the city of Barcelona. <br>


#### Privacy Policy 

In compliance with the obligations established in Organic Law 3/2018, of December 5, Protection of Personal Data and guarantee of digital rights and the General Regulation of Data Protection (RGPD), which is the Regulation (EU) 2016/679 of the Parliament and of the Council, of April 27, 2016, it is reported that the person in charge of the treatment of your personal data is Institut d'Investigació Biomètrica de Bellvitge (IDIBELL)CIF G58863317, with address Av/ Gran via de L'Hospitalet, nº 199-203, 08908, Hospitalet de Llobregat, his personal data included in the software (access IP), in order to manage its associated services.  <br>

The application only uses its session cookies for technical purposes. <br>

You will be able to request more information on the processing of data performed and exercise access rights, rectification, deletion, opposition, portability, and limitation of the process by sending a request through the email address dataprotection@idibell.cat. You will also be able to contact the data protection delegate of the IDIBELL via dpd@ticsalutsocial.cat. <br>  

Furthermore, if you consider that you have not obtained satisfaction in the exercise of your rights, you will be able to file a complaint before the Catalan Data Protection Agency, and you will be able to find more information on the processing of your data from the IDIBELL to the following Privacy Policy. <br>



