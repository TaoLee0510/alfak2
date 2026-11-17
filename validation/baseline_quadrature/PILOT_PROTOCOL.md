# Same-data deployment-prior development pilot

Use all twelve pre-existing finite-population process-development count pairs
(six landscapes, two count samples per landscape) from accuracy round 2.
They are revealed development data, not new validation data. The 27-node
domain, native observed-support distance, truth and fixed t0 reference remain
unchanged. No family label or truth enters the estimator.

Compare frozen V1 against the one unchanged baseline-integration algorithm.
Unlike the mechanism references' fixed SD-2 logit prior, this pilot preserves
the deployment initial-logit priors exactly: observed SD 5, borrowed mean -6,
borrowed SD 1.5 and distance penalty .75. A transparent fit bridge overrides
the mechanism fixture's SD-2 defaults, checks every conditional fit's returned
eta-prior vectors against frozen V1, and first replays V1's complete means and
covariance. The fitness prior and all quadrature criteria remain unchanged.

At most two fit processes run concurrently. Every positive quadrature node
must be trusted. Failed/nonconverged cases stay visible and cannot be dropped
from an aggregate. Both models use 40,000 joint distribution draws, with common
random numbers, for 50/80/95% reference-contrast intervals. The candidate samples
the actual component mixture, never a moment-matched log-Gaussian.

Report point RMSE, globally centered shape error, edge gradients and raw proper
interval scores at each native distance and by generating family. Average
sampling replicates within landscape, landscapes within family, then families
equally. No calibration is fitted here. Improvement in numerical posterior
agreement does not by itself establish improved generating-truth estimation.
The existing all-distance/family shape and uncertainty stage gates remain;
an inconclusive/failed pilot does not authorize opening reserved tests.
