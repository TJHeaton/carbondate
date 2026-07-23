# Outputs code suitable for running in OxCal from a series of radiocarbon determinations

Outputs code suitable for running in OxCal from a series of radiocarbon
determinations that can be given as either \\{}^{14}\\C age or
F\\{}^{14}\\C.

## Usage

``` r
GenerateOxcalCode(
  model_name,
  rc_determinations,
  rc_sigmas,
  rc_names = NULL,
  F14C_inputs = FALSE,
  outfile_path = NULL
)
```

## Arguments

- model_name:

  The name given to the model in the OxCal code.

- rc_determinations:

  A vector of observed radiocarbon determinations. Can be provided
  either as \\{}^{14}\\C ages (in \\{}^{14}\\C yr BP) or as
  F\\{}^{14}\\C concentrations.

- rc_sigmas:

  A vector of the (1-sigma) measurement uncertainties for the
  radiocarbon determinations. Must be the same length as
  `rc_determinations` and given in the same units.

- rc_names:

  Optional. The name of each data point - if given it must be the same
  length of rc_determinations.

- F14C_inputs:

  `TRUE` if the provided `rc_determinations` are F\\{}^{14}\\C
  concentrations and `FALSE` if they are radiocarbon ages. Defaults to
  `FALSE`.

- outfile_path:

  Optional. If given the OxCal code will be output to the file at the
  path given, otherwise it will be output to the terminal.

## Value

None

## Examples

``` r
# Generate names automatically and outputs to the screen for 14C ages
GenerateOxcalCode("My_data", c(1123, 1128, 1135), c(32, 24, 25))
#>  Plot()
#>  {
#>   NP_Model("My_data")
#>   {
#>   R_Date("1",1123,32);
#>   R_Date("2",1128,24);
#>   R_Date("3",1135,25);
#>   };
#>  };

# Provide name automatically and outputs to the screen for F14C concentrations
GenerateOxcalCode(
  "My_data",
  c(0.832, 0.850, 0.846),
  c(0.004, 0.003, 0.009),
  c("P-1", "P-2", "P-3"),
  F14C_inputs=TRUE)
#>  Plot()
#>  {
#>   NP_Model("My_data")
#>   {
#>   R_F14C("P-1",0.832,0.004);
#>   R_F14C("P-2",0.85,0.003);
#>   R_F14C("P-3",0.846,0.009);
#>   };
#>  };
```
