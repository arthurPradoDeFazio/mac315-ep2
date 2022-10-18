0;
function matrix = B_barinv(Binv, m, u, l) # monta a matriz [B^-1 | u] e faz com que u se transforme no l-ésimo vetor canonico, o que está à esquerda é B barra inversa
  matrix = [Binv u];
  n_rows = m;
  n_columns = m + 1;
  
  for i = 1:n_rows
    if i == l
      matrix(l, :) = (1/matrix(l, n_columns)) * matrix(l, :);
    else
      matrix(i, :) = matrix(i, :) - (matrix(i, n_columns) / matrix(l, n_columns)) * matrix(l, :);
    endif
  endfor
  
  matrix(:, [n_columns]) = [];
endfunction

function [theta_star, l] = theta(x, bind, m, u) # calcula o theta para andar na direcão básica e o índice em que ele é atingido
  theta_star = Inf;
  l = -1;
  
  for i = 1:m
    if u(i) > 0
      new_theta = x(bind(i)) / u(i);
      if new_theta < theta_star # note que escolhemos o primeiro indice l que realiza o mínimo
        theta_star = new_theta;
        l = i;
      endif
    endif
  endfor
endfunction

function j_index = entering_column_index(A, c, bind, Binv, m, n) # talve esteja errada
  bind_sorted = sort(bind);
  looking_at = 1;
  p_transpose = c(bind)' * Binv;
  j_index = -1;
  
  for j = 1:n
    # se j nao pertence ao vetor bind
    if looking_at > m || j != bind_sorted(looking_at) 
      if c(j) - p_transpose * A(:, j) < 0 # escolhemos o primeiro custo negativo
        j_index = j;
        return;
      endif
    else
      looking_at += 1;
    endif
  endfor
  
endfunction

function y = new_basic_solution(x, bind, m, n, u, theta, new_basic_index)
  y = zeros(n, 1);
  y(new_basic_index) = theta;
  for i = 1:m
    y(bind(i)) = x(bind(i)) - theta * u(i);
  endfor
endfunction

function [ind v] = simplex(A, b, c, m, n, x, bind, Binv)
  while 1
    new_basic_index = entering_column_index(A, c, bind, Binv, m, n);
    if new_basic_index == -1
      ind = 0;
      v = x;
      return;
    endif
    
    u = Binv * A(:, new_basic_index);
    [theta, exiting_index] = theta_star(x, bind, m, u)
    if theta == Inf
      ind = 0;
      v = -u;
      return;
    endif
    
    bind(exiting_index) = new_basic_index;
    x = new_basic_solution(x, bind, m, n, u, theta, new_basic_index);
    Binv = B_barinv(Binv, m, u, l);
  endwhile
endfunction
