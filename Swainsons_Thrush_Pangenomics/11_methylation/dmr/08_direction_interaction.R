# Is the hypo- vs hypermethylated asymmetry in TE overlap stronger in hybrid comparisons than between the
# pure subspecies? Binomial GLM of TE overlap on DMR direction, comparison type (pure / hybrid) and their
# interaction, with the interaction tested by a likelihood ratio test against the additive model.

library(data.table)

DIR <- "methylation/dmr/te_enrichment"
d <- fread(file.path(DIR, "dmr_direction_fisher.tsv"))[annotation %in% c("EDTA", "RepeatModeler", "GraffiTE")]

long <- rbind(d[, .(comparison, annotation, direction = "hypo", hit = hypo_hit, n = hypo_n)],
              d[, .(comparison, annotation, direction = "hyper", hit = hyper_hit, n = hyper_n)])
long[, `:=`(direction = factor(direction, c("hyper", "hypo")),
            type = factor(fifelse(comparison == "coastal_vs_inland", "pure", "hybrid"), c("pure", "hybrid")))]

res <- rbindlist(lapply(split(long, by = "annotation"), function(x) {
    full <- glm(cbind(hit, n - hit) ~ direction * type, family = binomial, data = x)
    add  <- glm(cbind(hit, n - hit) ~ direction + type, family = binomial, data = x)
    b <- coef(full)
    data.table(annotation = x$annotation[1],
               or_pure = exp(b[["directionhypo"]]),
               or_hybrid = exp(b[["directionhypo"]] + b[["directionhypo:typehybrid"]]),
               lrt_p = anova(add, full, test = "LRT")[2, "Pr(>Chi)"])
}))
fwrite(res, file.path(DIR, "direction_interaction.tsv"), sep = "\t")
