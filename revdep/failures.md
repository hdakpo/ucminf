# mets

<details>

* Version: 1.3.2
* GitHub: https://github.com/kkholst/mets
* Source code: https://github.com/cran/mets
* Date/Publication: 2023-01-17 09:40:07 UTC
* Number of recursive dependencies: 95

Run `revdepcheck::revdep_details(, "mets")` for more info

</details>

## In both

*   R CMD check timed out
    

*   checking C++ specification ... NOTE
    ```
      Specified C++11: please drop specification unless essential
    ```

*   checking installed package size ... NOTE
    ```
      installed size is  8.4Mb
      sub-directories of 1Mb or more:
        R      1.7Mb
        doc    2.6Mb
        libs   2.3Mb
    ```

# nlmixr2est

<details>

* Version: 2.1.5
* GitHub: https://github.com/nlmixr2/nlmixr2est
* Source code: https://github.com/cran/nlmixr2est
* Date/Publication: 2023-04-22 19:50:02 UTC
* Number of recursive dependencies: 198

Run `revdepcheck::revdep_details(, "nlmixr2est")` for more info

</details>

## In both

*   checking tests ...
    ```
      Running 'testthat.R'
     ERROR
    Running the tests in 'tests/testthat.R' failed.
    Last 13 lines of output:
      
      ══ Failed tests ════════════════════════════════════════════════════════════════
      ── Error ('test-00-reload-ll.R:30:3'): between session saem ll works ───────────
      Error in `gzfile(file, "rb")`: cannot open the connection
      Backtrace:
          ▆
       1. ├─withr::with_tempdir(...) at test-00-reload-ll.R:30:2
       2. │ └─withr::with_dir(tmp, code)
       3. │   └─base::force(code)
       4. └─base::readRDS("fit.rds") at test-00-reload-ll.R:36:4
       5.   └─base::gzfile(file, "rb")
      
      [ FAIL 1 | WARN 1 | SKIP 0 | PASS 157 ]
      Error: Test failures
      Execution halted
    ```

*   R CMD check timed out
    

