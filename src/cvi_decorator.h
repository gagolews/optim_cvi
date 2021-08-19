#ifndef __CVI_DECORATOR_H
#define __CVI_DECORATOR_H

#include "cvi.h"
#include "cvi_external.h"


class ClusterValidityIndexDecorator
{

public:
    virtual void modify_with_weight(size_t i, uint8_t j, size_t w, const std::unordered_map<size_t, size_t>& new2old) = 0;
    virtual void set_labels_with_weights(const std::vector<uint8_t>& _L, const std::vector<size_t>& _weights, const std::unordered_map<size_t, size_t>& new2old) = 0;
    virtual void set_labels(const std::vector<uint8_t>& _L) = 0;
    virtual uint8_t get_label_translated(size_t i, const std::unordered_map<size_t, size_t>& new2old) = 0;
    virtual std::vector<uint8_t> get_labels_translated(size_t n, const std::unordered_map<size_t, size_t>& new2old) = 0;
    virtual FLOAT_T compute() = 0;
};


class InternalClusterValidityIndexDecorator : public ClusterValidityIndexDecorator
{
    ClusterValidityIndex* index;
public:
    InternalClusterValidityIndexDecorator(ClusterValidityIndex* index)
        : index(index)
    {

    }

    virtual void modify_with_weight(size_t i, uint8_t j, size_t w, const std::unordered_map<size_t, size_t>& new2old)
    {
        index->modify(i, j);
    }

    virtual void set_labels_with_weights(const std::vector<uint8_t>& _L, const std::vector<size_t>& _weights, const std::unordered_map<size_t, size_t>& new2old)
    {
        index->set_labels(_L);
    }

    virtual void set_labels(const std::vector<uint8_t>& _L)
    {
        index->set_labels(_L);
    }

    virtual uint8_t get_label_translated(size_t i, const std::unordered_map<size_t, size_t>& new2old)
    {
//        Rcpp::Rcout << "w get_label_translated" << std::endl;
        return index->get_label(i);
    }

    virtual std::vector<uint8_t> get_labels_translated(size_t n, const std::unordered_map<size_t, size_t>& new2old)
    {
        return index->get_labels();
    }

    virtual FLOAT_T compute()
    {
        return index->compute();
    }
};

class ExternalClusterValidityIndexDecorator : public ClusterValidityIndexDecorator
{
    ExternalClusterValidityIndex* index;
public:
    ExternalClusterValidityIndexDecorator(ExternalClusterValidityIndex* index)
        : index(index)
    {

    }

    virtual void modify_with_weight(size_t i, uint8_t j, size_t w, const std::unordered_map<size_t, size_t>& new2old)
    {
        index->modify_with_weight(i, j, w, new2old);
    }

    virtual void set_labels_with_weights(const std::vector<uint8_t>& _L, const std::vector<size_t>& _weights, const std::unordered_map<size_t, size_t>& new2old)
    {
        index->set_labels_with_weights(_L, _weights, new2old);
    }

    virtual void set_labels(const std::vector<uint8_t>& _L)
    {
        index->set_labels(_L);
    }

    virtual uint8_t get_label_translated(size_t i, const std::unordered_map<size_t, size_t>& new2old)
    {
        return index->get_label(new2old.at(i));
    }

    virtual std::vector<uint8_t> get_labels_translated(size_t n, const std::unordered_map<size_t, size_t>& new2old)
    {
        std::vector<uint8_t> result(n);
        for(size_t i = 0; i < n; ++i) {
            result[i] = index->get_label(new2old.at(i));
        }
        return result;
    }

    virtual FLOAT_T compute()
    {
        return index->compute();
    }
};

#endif
