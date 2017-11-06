function c = load_project_const()
  % given constants
  c = struct();

  c.v_dd = 1.2;
  c.l_min = 90e-9;
  c.w_min = 180e-9;

  c.snr_min_db = 55;

  c.p_meta = 10e-7;

  c.bits = 10;
  
  c.min_cap = 2.5e-15;
  c.c_in_sar_logic = 50e-15;
  c.v_ic_min = 0.4;
  c.v_ic_max = 0.8;
  
  c.bias_current = 100e-6;
  
  % derived constants
  
  c.n_codes = 2 ^ c.bits;
  
end
