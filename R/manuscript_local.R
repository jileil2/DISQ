# Decorrelated roots with each worker's fitted SQBiT nuisance held fixed.
.local_dsqbit_need <- function(ok,message) if(!isTRUE(ok))stop(message,call.=FALSE)
.local_dsqbit_fail <- function(message,worker,diagnostics=NULL) {
  stop(structure(list(message=paste0("Local dSQBiT worker ",worker,": ",message),
    call=NULL,worker=worker,diagnostics=diagnostics),
    class=c("local_dsqbit_error","error","condition")))
}

.disq_local_roots <- function(Ys,Cs,Zs,nuisance,eta_local,phi_local,tau,
                                     kernel="norm",tolerance=1e-11,max_iter=200L,
                                     raw_residuals=NULL) {
  kernel<-match.arg(kernel,c("norm","unif"))
  .local_dsqbit_need(is.numeric(tau)&&length(tau)==1L&&is.finite(tau)&&tau>0&&tau<1,"Invalid tau")
  .local_dsqbit_need(is.numeric(tolerance)&&length(tolerance)==1L&&is.finite(tolerance)&&tolerance>0&&
    is.numeric(max_iter)&&length(max_iter)==1L&&is.finite(max_iter)&&max_iter>=1&&max_iter==floor(max_iter),
    "Invalid local score solve controls")
  K<-length(Ys);eta_local<-as.matrix(eta_local);q<-ncol(eta_local)
  .local_dsqbit_need(K>0&&length(Cs)==K&&length(Zs)==K&&length(nuisance)==K&&
    nrow(eta_local)==K&&q>0&&all(is.finite(eta_local))&&length(phi_local)==K&&
    all(is.finite(phi_local)&phi_local>0),"Saved local input dimensions or bandwidths are invalid")
  if(!is.null(raw_residuals)).local_dsqbit_need(length(raw_residuals)==K,"Saved raw residual roster differs")
  eta_star<-matrix(NA_real_,K,q,dimnames=dimnames(eta_local))
  diagnostics<-vector("list",K);local_roots<-vector("list",K);projections<-vector("list",K)
  for(k in seq_len(K)) {
    completed<-function()if(k>1L)do.call(rbind,diagnostics[seq_len(k-1L)]) else NULL
    one<-tryCatch({
      Yk<-as.numeric(Ys[[k]]);Ck<-as.matrix(Cs[[k]]);Zk<-as.matrix(Zs[[k]])
      muk<-as.numeric(nuisance[[k]]);eta0<-as.numeric(eta_local[k,]);phi<-as.numeric(phi_local[k]);nk<-length(Yk)
      .local_dsqbit_need(nk>0&&nrow(Ck)==nk&&ncol(Ck)==q&&nrow(Zk)==nk&&ncol(Zk)>0&&
        length(muk)==nk&&all(is.finite(c(Yk,Ck,Zk,muk))),"Nonfinite or mismatched local observations/design")
      raw_residual<-Yk-as.numeric(Ck%*%eta0)-muk
      residual_difference<-NA_real_
      if(!is.null(raw_residuals)) {
        saved<-as.numeric(raw_residuals[[k]])
        .local_dsqbit_need(length(saved)==nk&&all(is.finite(saved)),"Invalid saved raw SQBiT residual")
        residual_difference<-max(abs(raw_residual-saved))
        .local_dsqbit_need(residual_difference<=1e-8*max(1,max(abs(saved))),
          "Reconstructed LOCAL nuisance/design does not reproduce the saved raw SQBiT residual")
      }
      # Exactly the existing SQBiT_ds construction: local raw nuisance and raw eta
      # determine the fixed full-spline density weights at the saved local phi.
      projection<-.disq_local_projection(Zk,Ck,raw_residual,phi,kernel=kernel)
      .local_dsqbit_need(projection$rank==ncol(Zk)&&
        is.finite(projection$orthogonality_error)&&projection$orthogonality_error<=1e-8,
        "The local full-spline projection did not pass rank/normal-equation checks")
      fit<-.disq_solve_score(eta0,Yk,Ck,muk,projection$residualized,phi,tau,
        max.iter=as.integer(max_iter),tolerance=tolerance,kernel=kernel)
      .local_dsqbit_need(isTRUE(fit$converged)&&length(fit$eta)==q&&all(is.finite(fit$eta))&&
        is.finite(fit$score_norm)&&fit$score_norm<=tolerance,"Local decorrelated score solve did not converge")
      fitted_residual<-Yk-muk-as.numeric(Ck%*%fit$eta)
      # Preserve the nonsymmetric D'C orientation; do not replace C with D.
      density<-if(kernel=="norm")dnorm(fitted_residual/phi)/phi else dunif(fitted_residual/phi,-1,1)/phi
      J<-crossprod(projection$residualized,Ck*density)/nk
      source_J<-.disq_jacobian(fit$eta,Yk,Ck,muk,projection$residualized,phi,tau,kernel=kernel)
      J_difference<-max(abs(J-source_J))
      J_rcond<-rcond(J)
      .local_dsqbit_need(all(is.finite(J))&&is.finite(J_rcond)&&J_rcond>0&&
        J_difference<=1e-12*max(1,max(abs(J))),"Local D'C Jacobian is invalid or differs from the accepted score derivative")
      distribution<-if(kernel=="norm")pnorm(-fitted_residual/phi) else punif(-fitted_residual/phi,-1,1)
      score<-as.numeric(colMeans(projection$residualized*(distribution-tau)))
      .local_dsqbit_need(max(abs(score))<=tolerance&&max(abs(score-fit$score))<=1e-12,
        "Independent local score readback differs from the solved equation")
      row<-data.frame(worker=k,n=nk,nuisance_dimension=ncol(Zk),phi=phi,
        converged=TRUE,iterations=fit$iterations,score_norm=max(abs(score)),tolerance=tolerance,
        max_iter=as.integer(max_iter),projection_rank=projection$rank,
        projection_condition=projection$condition,orthogonality_error=projection$orthogonality_error,
        raw_residual_reproduction_max_difference=residual_difference,
        jacobian_rcond=J_rcond,jacobian_asymmetry=max(abs(J-t(J))),
        jacobian_source_max_difference=J_difference,
        eta_change_l2=sqrt(sum((fit$eta-eta0)^2)),root_method=fit$method)
      # Compact outputs: no observation-level designs or residualized arrays are retained.
      list(eta=as.numeric(fit$eta),diagnostics=row,
        root=c(fit,list(J=J)),projection=projection[c("P","rank","condition","orthogonality_error")])
    },error=function(error).local_dsqbit_fail(conditionMessage(error),k,completed()))
    eta_star[k,]<-one$eta;diagnostics[[k]]<-one$diagnostics
    local_roots[[k]]<-one$root;projections[[k]]<-one$projection
  }
  list(initializer_type="unweighted_mean_local_dSQBiT",eta_star=eta_star,
    eta_mean=as.numeric(colMeans(eta_star)),raw_eta=eta_local,raw_mean=as.numeric(colMeans(eta_local)),
    diagnostics=do.call(rbind,diagnostics),local_roots=local_roots,projections=projections,
    kernel=kernel,tau=tau,tolerance=tolerance,max_iter=as.integer(max_iter),
    nuisance_policy="Each original LOCAL raw SQBiT nuisance Z_k*saved gamma_k is held fixed",
    projection_weights="Saved local phi and original raw SQBiT eta/nuisance residuals",
    original_initial_fit_modified=FALSE)
}
