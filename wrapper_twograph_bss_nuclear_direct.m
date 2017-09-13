function wrapper_twograph_bss_nuclear_direct

verbose_twograph_bss_nuclear_direct = true;
verbose_self = true;

[truth, model, y] = twograph_bss_gen_problem;
[Z1_hat, Z2_hat] = twograph_bss_nuclear_direct(y, model.A1, model.G(1).V, ...
                                               model.A2, model.G(2).V, ...
                                               verbose_twograph_bss_nuclear_direct);
if verbose_self
  twograph_bss_print_summary(Z1_hat, Z2_hat, truth, model, y);
end

end
