# mets

<details>

* Version: 1.3.4
* GitHub: https://github.com/kkholst/mets
* Source code: https://github.com/cran/mets
* Date/Publication: 2024-02-16 19:10:07 UTC
* Number of recursive dependencies: 95

Run `revdepcheck::revdep_details(, "mets")` for more info

</details>

## Newly broken

*   R CMD check timed out
    

## In both

*   checking running R code from vignettes ...
    ```
      'basic-dutils.Rmd' using 'UTF-8'... OK
      'binomial-family.Rmd' using 'UTF-8'... OK
      'binomial-twin.Rmd' using 'UTF-8'... OK
      'binreg-TRS.rmd' using 'UTF-8'... OK
      'binreg-ate.Rmd' using 'UTF-8'... OK
      'binreg.Rmd' using 'UTF-8'... OK
      'cifreg.Rmd' using 'UTF-8'... OK
      'glm-utility.Rmd' using 'UTF-8'... OK
      'haplo-discrete-ttp.Rmd' using 'UTF-8'... OK
      'interval-discrete-survival.Rmd' using 'UTF-8'... OK
    ...
    > margph <- phreg(Surv(time, cancer) ~ strata(country) + 
    +     cluster(id), data = prt)
    
    > plot(margph)
    
    > summary(fitco1)
    
      When sourcing 'time-to-event-family-studies-arev.R':
    Erreur : objet 'fitco1' introuvable
    Exécution arrêtée
    ```

*   checking installed package size ... NOTE
    ```
      installed size is  8.7Mb
      sub-directories of 1Mb or more:
        R      1.7Mb
        doc    2.8Mb
        libs   2.3Mb
    ```

