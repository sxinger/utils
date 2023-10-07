************************************
## Markdown page setup chunk

```r
knitr::opts_chunk$set(echo = FALSE,    # hide codes
                      message = FALSE, # hide message
                      warning = FALSE, # hide warning 
                      fig.height=10,
                      fig.width = 6
                      )
```


************************************
## Footer

*Copyright (c) 2023 The Curators of University of Missouri* <br/>
*Share and Enjoy according to the terms of the MIT Open Source License* <br/>
*development conventions and git collaborative workflow: https://github.com/RWD2E* <br/>
*Repo structure style follows [GitHub Repository Structure Best Practices](https://soulaimanghanem.medium.com/github-repository-structure-best-practices-248e6effc405) (Ghanem, 2021)*


****************************************
## Headers

### SQL script header
/*
# Copyright (c) 2021-2025 University of Missouri                   
# Author: Xing Song, xsm7f@umsystem.edu                            
# File: <script_name>.sql
# Description: short description of the script
*/

### R script headers
#################################################################
# Copyright (c) 2021-2025 University of Missouri                   
# Author: Xing Song, xsm7f@umsystem.edu                            
# File: extract.R
# Description: pre-analytic data from snowflake to R
#################################################################


### Snowflake stored procedures templates
create or replace procedure foo(
    ...
    DRY_RUN boolean,
    DRY_RUN_AT string
)
returns variant
language javascript
as
$$
/**
 * @param {array} ...: other parameters
 * @param {boolean} DRY_RUN: dry run indicator. If true, only sql script will be created and stored in dev.sp_out table
 * @param {boolean} DRY_RUN_AT: A temporary location to store the generated sql query for debugging purpose. 
                                When DRY_RUN = True, provide absolute path to the table; when DRY_RUN = False, provide NULL 
**/
if (DRY_RUN) {
    var log_stmt = snowflake.createStatement({
        sqlText: `CREATE OR REPLACE TEMPORARY TABLE `+ DRY_RUN_AT +`(QRY VARCHAR);`});
    log_stmt.execute(); 
}

... main body ... 

if (DRY_RUN) {
        // preview of the generated dynamic SQL scripts - comment it out when perform actual execution
        var log_stmt = snowflake.createStatement({
                        sqlText: `INSERT INTO `+ DRY_RUN_AT +` (qry) values (:1), (:2);`,
                        binds: [sqlstmt1,sqlstmt2]});
        log_stmt.execute(); 
    } else {
        // run dynamic dml query
        var run_sqlstmt1 = snowflake.createStatement({sqlText: sqlstmt1}); run_sqlstmt1.execute();
        var run_sqlstmt2 = snowflake.createStatement({sqlText: sqlstmt2}); run_sqlstmt2.execute();
        var commit_txn = snowflake.createStatement({sqlText: `commit;`}); commit_txn.execute();
    }

$$
;