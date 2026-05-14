// src/marley_hepmc3.cc
//
// Amalgamated built-in HepMC3 implementation for MARLEY.
//
// This translation unit is compiled when no system HepMC3 installation is
// available (MARLEY_FOUND_HEPMC3 is not defined by the build system). When
// a system installation is used, the preprocessor guard below renders this
// file a no-op; it is still compiled and linked into libHepMC3 but
// contributes no symbols.
//
// Provenance
// ----------
// Derived from the HepMC3 project: https://gitlab.cern.ch/hepmc/HepMC3
// Authors: the HepMC3 collaboration.
// Redistributed under the GNU Lesser General Public License, version 2.1
// or later. Copyright (C) 2014-2023 the HepMC Collaboration.
//
// Only the subset of source files required for MARLEY's output support is
// included (in compilation order):
//   hepmc3/Setup.cc
//   hepmc3/GenRunInfo.cc
//   hepmc3/GenParticle.cc
//   hepmc3/GenVertex.cc
//   hepmc3/GenEvent.cc
//   hepmc3/Print.cc
//   hepmc3/WriterAscii.cc
//   hepmc3/ReaderAscii.cc
// The bundled header files live in include/builtin/HepMC3/ (see README there).

#ifndef MARLEY_FOUND_HEPMC3

// ============================================================
// hepmc3/Setup.cc
// Copyright (C) 2014-2023 The HepMC collaboration
// ============================================================
// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2023 The HepMC collaboration (see AUTHORS for details)
//
/**
 *  @file Setup.cc
 *  @brief Implementation of \b Setup class
 *
 */
#include "HepMC3/Setup.h"

namespace HepMC3 {

const unsigned int Setup::DEFAULT_DOUBLE_ALMOST_EQUAL_MAXULPS = 10;
const double       Setup::DOUBLE_EPSILON                      = 10e-20;
bool Setup::print_errors()                      { return m_is_printing_errors;    }
void Setup::set_print_errors(const bool flag)   { m_is_printing_errors   = flag;  }
bool Setup::print_warnings()                    { return m_is_printing_warnings;  }
void Setup::set_print_warnings(const bool flag) { m_is_printing_warnings = flag;  }
int  Setup::debug_level()                       { return m_debug_level;           }
void Setup::set_debug_level(const int level)    { m_debug_level          = level; }
int  Setup::errors_level()                      { return m_errors_level;           }
void Setup::set_errors_level(const int level)   { m_errors_level          = level;  }
int  Setup::warnings_level()                    { return m_warnings_level;          }
void Setup::set_warnings_level(const int level) { m_warnings_level          = level;}
bool Setup::m_is_printing_errors    = true;
bool Setup::m_is_printing_warnings  = true;
int  Setup::m_debug_level           = 5;
int  Setup::m_errors_level           = 1000;
int  Setup::m_warnings_level           = 750;

} // namespace HepMC3

// ============================================================
// hepmc3/GenRunInfo.cc
// Copyright (C) 2014-2023 The HepMC collaboration
// ============================================================
// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2023 The HepMC collaboration (see AUTHORS for details)
//
/**
 *  @file GenRunInfo.cc
 *  @brief Implementation of \b class GenRunInfo
 *
 */
#include <sstream>

#include "HepMC3/Data/GenRunInfoData.h"
#include "HepMC3/GenRunInfo.h"


namespace HepMC3 {


void GenRunInfo::set_weight_names(const std::vector<std::string> & names) {
    m_weight_indices.clear();
    m_weight_names = names;
    for ( int i = 0, N = names.size(); i < N; ++i ) {
        std::string name = names[i];
        if ( name.empty() ) {
            std::ostringstream oss;
            oss << i;
            name = oss.str();
            m_weight_names[i] = name;
        }
        if ( has_weight(name) ) {
            throw std::logic_error("GenRunInfo::set_weight_names: "
                                   "Duplicate weight name '" + name +
                                   "' found.");
        }
        m_weight_indices[name] = i;
    }
}

std::string GenRunInfo::attribute_as_string(const std::string &name) const {
    std::lock_guard<std::recursive_mutex> lock(m_lock_attributes);
    auto i = m_attributes.find(name);
    if ( i == m_attributes.end() ) return {};

    if ( !i->second ) return {};

    std::string ret;
    i->second->to_string(ret);

    return ret;
}

void GenRunInfo::write_data(GenRunInfoData& data) const {
    // Weight names
    data.weight_names = this->weight_names();

    // Attributes
    using att_val_t = std::map<std::string, std::shared_ptr<Attribute>>::value_type;

    for (const att_val_t& vt: m_attributes) {
        std::string att;
        vt.second->to_string(att);

        data.attribute_name.  emplace_back(vt.first);
        data.attribute_string.emplace_back(att);
    }

    // Tools
    for ( const ToolInfo &tool: this->tools() ) {
        data.tool_name.       emplace_back(tool.name);
        data.tool_version.    emplace_back(tool.version);
        data.tool_description.emplace_back(tool.description);
    }
}


std::vector<std::string> GenRunInfo::attribute_names() const {
    std::vector<std::string> results;
    results.reserve(m_attributes.size());
    for (const auto& vt1: m_attributes) {
        results.emplace_back(vt1.first);
    }
    return results;
}

void GenRunInfo::read_data(const GenRunInfoData& data) {
    // Weight names
    set_weight_names(data.weight_names);

    // Attributes
    for (unsigned int i = 0; i < data.attribute_name.size(); ++i) {
        add_attribute(data.attribute_name[i],
                      std::make_shared<StringAttribute>(data.attribute_string[i]));
    }

    // Tools
    for (unsigned int i = 0; i < data.tool_name.size(); ++i) {
        ToolInfo ti;
        ti.name        = data.tool_name[i];
        ti.version     = data.tool_version[i];
        ti.description = data.tool_description[i];

        this->tools().emplace_back(ti);
    }
}

GenRunInfo::GenRunInfo(const GenRunInfo& r)
{
    if (this != &r)
    {
        std::lock(m_lock_attributes, r.m_lock_attributes);
        std::lock_guard<std::recursive_mutex> lhs_lk(m_lock_attributes, std::adopt_lock);
        std::lock_guard<std::recursive_mutex> rhs_lk(r.m_lock_attributes, std::adopt_lock);
        GenRunInfoData tdata;
        r.write_data(tdata);
        read_data(tdata);
    }
}
GenRunInfo& GenRunInfo::operator=(const GenRunInfo& r)
{
    if (this != &r)
    {
        std::lock(m_lock_attributes, r.m_lock_attributes);
        std::lock_guard<std::recursive_mutex> lhs_lk(m_lock_attributes, std::adopt_lock);
        std::lock_guard<std::recursive_mutex> rhs_lk(r.m_lock_attributes, std::adopt_lock);
        GenRunInfoData tdata;
        r.write_data(tdata);
        read_data(tdata);
    }
    return *this;
}

} // namespace HepMC3

// ============================================================
// hepmc3/GenParticle.cc
// Copyright (C) 2014-2023 The HepMC collaboration
// ============================================================
// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2023 The HepMC collaboration (see AUTHORS for details)
//
/**
 *  @file GenParticle.cc
 *  @brief Implementation of \b class GenParticle
 *
 */
#include "HepMC3/Attribute.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Setup.h"


namespace HepMC3 {

GenParticle::GenParticle(const FourVector &mom, int pidin, int stat):
    m_event(nullptr),
    m_id(0) {
    m_data.pid               = pidin;
    m_data.momentum          = mom;
    m_data.status            = stat;
    m_data.is_mass_set       = false;
    m_data.mass              = 0.0;
}

GenParticle::GenParticle(const GenParticleData &dat):
    m_event(nullptr),
    m_id(0),
    m_data(dat) {
}

double GenParticle::generated_mass() const {
    return m_data.is_mass_set ? m_data.mass : m_data.momentum.m();
}

void GenParticle::set_pid(int pidin) {
    m_data.pid = pidin;
}

void GenParticle::set_status(int stat) {
    m_data.status = stat;
}

void GenParticle::set_momentum(const FourVector& mom) {
    m_data.momentum = mom;
}

void GenParticle::set_generated_mass(double m) {
    m_data.mass        = m;
    m_data.is_mass_set = true;
}

void GenParticle::unset_generated_mass() {
    m_data.mass        = 0.;
    m_data.is_mass_set = false;
}

GenVertexPtr GenParticle::production_vertex() {
    return m_production_vertex.lock();
}

ConstGenVertexPtr GenParticle::production_vertex() const {
    return std::const_pointer_cast<const GenVertex>(m_production_vertex.lock());
}

GenVertexPtr GenParticle::end_vertex() {
    return m_end_vertex.lock();
}

ConstGenVertexPtr GenParticle::end_vertex() const {
    return std::const_pointer_cast<const GenVertex>(m_end_vertex.lock());
}

std::vector<GenParticlePtr> GenParticle::parents() {
    return (m_production_vertex.expired())? std::vector<GenParticlePtr>() : production_vertex()->particles_in();
}

std::vector<ConstGenParticlePtr> GenParticle::parents() const {
    return (m_production_vertex.expired()) ? std::vector<ConstGenParticlePtr>() : production_vertex()->particles_in();
}

std::vector<GenParticlePtr> GenParticle::children() {
    return (m_end_vertex.expired())? std::vector<GenParticlePtr>() : end_vertex()->particles_out();
}

std::vector<ConstGenParticlePtr> GenParticle::children() const {
    return (m_end_vertex.expired()) ? std::vector<ConstGenParticlePtr>() : end_vertex()->particles_out();
}

bool GenParticle::add_attribute(const std::string& name, std::shared_ptr<Attribute> att) {
    if ( !parent_event() ) return false;
    parent_event()->add_attribute(name, att, id());
    return true;
}

std::vector<std::string> GenParticle::attribute_names() const {
    if ( parent_event() ) return parent_event()->attribute_names(id());
    return {};
}

void GenParticle::remove_attribute(const std::string& name) {
    if ( parent_event() ) parent_event()->remove_attribute(name, id());
}

std::string GenParticle::attribute_as_string(const std::string& name) const {
    return parent_event() ? parent_event()->attribute_as_string(name, id()) : std::string();
}

} // namespace HepMC3

// ============================================================
// hepmc3/GenVertex.cc
// Copyright (C) 2014-2023 The HepMC collaboration
// ============================================================
// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2023 The HepMC collaboration (see AUTHORS for details)
//
/**
 *  @file GenVertex.cc
 *  @brief Implementation of \b class GenVertex
 *
 */
#include <algorithm> // std::remove

#include "HepMC3/Attribute.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Setup.h"

namespace HepMC3 {


GenVertex::GenVertex(const FourVector& pos):
    m_event(nullptr),
    m_id(0) {
    m_data.status   = 0;
    m_data.position = pos;
}

GenVertex::GenVertex(const GenVertexData &dat):
    m_event(nullptr),
    m_id(0),
    m_data(dat) {
}


void GenVertex::add_particle_in(GenParticlePtr p) {
    if (!p) return;

    // Avoid duplicates
    if (std::find(particles_in().begin(), particles_in().end(), p) != particles_in().end()) return;

    m_particles_in.emplace_back(p);

    if ( p->end_vertex() ) p->end_vertex()->remove_particle_in(p);

    p->m_end_vertex = shared_from_this();

    if (m_event) m_event->add_particle(p);
}


void GenVertex::add_particle_out(GenParticlePtr p) {
    if (!p) return;

    // Avoid duplicates
    if (std::find(particles_out().begin(), particles_out().end(), p) != particles_out().end()) return;

    m_particles_out.emplace_back(p);

    if ( p->production_vertex() ) p->production_vertex()->remove_particle_out(p);

    p->m_production_vertex = shared_from_this();

    if (m_event) m_event->add_particle(p);
}

void GenVertex::remove_particle_in(GenParticlePtr p) {
    if (!p) return;
    if (std::find(m_particles_in.begin(), m_particles_in.end(), p) == m_particles_in.end()) return;
    p->m_end_vertex.reset();
    m_particles_in.erase(std::remove(m_particles_in.begin(), m_particles_in.end(), p), m_particles_in.end());
}


void GenVertex::remove_particle_out(GenParticlePtr p) {
    if (!p) return;
    if (std::find(m_particles_out.begin(), m_particles_out.end(), p) == m_particles_out.end()) return;
    p->m_production_vertex.reset();
    m_particles_out.erase(std::remove(m_particles_out.begin(), m_particles_out.end(), p), m_particles_out.end());
}

const std::vector<ConstGenParticlePtr>& GenVertex::particles_in()const {
    return *(reinterpret_cast<const std::vector<ConstGenParticlePtr>*>(&m_particles_in));
}

const std::vector<ConstGenParticlePtr>& GenVertex::particles_out()const {
    return *(reinterpret_cast<const std::vector<ConstGenParticlePtr>*>(&m_particles_out));
}

const FourVector& GenVertex::position() const {
    if ( has_set_position() ) return m_data.position;

    // No position information - look at event and/or search ancestors
    if ( parent_event() )
    {
        std::shared_ptr<IntAttribute> cycles = parent_event()->attribute<IntAttribute>("cycles");
        //This could be a recussive call.  Try to prevent it.
        if (!cycles || cycles->value() == 0)
        {
            for (const auto& p: m_particles_in) {
                ConstGenVertexPtr v = p->production_vertex();
                if (v) return v->position();
            }
        }
        return parent_event()->event_pos();
    }
    return FourVector::ZERO_VECTOR();
}

void GenVertex::set_position(const FourVector& new_pos) {
    m_data.position = new_pos;
}

bool GenVertex::add_attribute(const std::string& name, std::shared_ptr<Attribute> att) {
    if ( !parent_event() ) return false;
    parent_event()->add_attribute(name, att, id());
    return true;
}

void GenVertex::remove_attribute(const std::string& name) {
    if ( parent_event() ) parent_event()->remove_attribute(name, id());
}

std::string GenVertex::attribute_as_string(const std::string& name) const {
    return parent_event() ? parent_event()->attribute_as_string(name, id()) : std::string();
}

std::vector<std::string> GenVertex::attribute_names() const {
    if ( parent_event() ) return parent_event()->attribute_names(id());

    return {};
}

} // namespace HepMC3

// ============================================================
// hepmc3/GenEvent.cc
// Copyright (C) 2014-2023 The HepMC collaboration
// ============================================================
// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2023 The HepMC collaboration (see AUTHORS for details)
//
/**
 *  @file GenEvent.cc
 *  @brief Implementation of \b class GenEvent
 *
 */
#include <algorithm> // sort
#include <deque>

#include "HepMC3/Data/GenEventData.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"


namespace HepMC3 {

GenEvent::GenEvent(Units::MomentumUnit mu,
                   Units::LengthUnit lu)
    : m_momentum_unit(mu), m_length_unit(lu), //m_weights(std::vector<double>(1, 1.0)),//Prevent from  different number of weights and names
      m_rootvertex(std::make_shared<GenVertex>()) {}


GenEvent::GenEvent(std::shared_ptr<GenRunInfo> run,
                   Units::MomentumUnit mu,
                   Units::LengthUnit lu)
    : m_momentum_unit(mu), m_length_unit(lu),  //m_weights(std::vector<double>(1, 1.0)),//Prevent from  different number of weights and names
      m_rootvertex(std::make_shared<GenVertex>()),
      m_run_info(run) {
    if ( run && !run->weight_names().empty() ) {
        m_weights = std::vector<double>(run->weight_names().size(), 1.0);
    }
}

const std::vector<ConstGenParticlePtr>& GenEvent::particles() const {
    return *(reinterpret_cast<const std::vector<ConstGenParticlePtr>*>(&m_particles));
}

const std::vector<ConstGenVertexPtr>& GenEvent::vertices() const {
    return *(reinterpret_cast<const std::vector<ConstGenVertexPtr>*>(&m_vertices));
}


void GenEvent::add_particle(GenParticlePtr p) {
    if ( !p || p->in_event() ) return;

    m_particles.emplace_back(p);

    p->m_event = this;
    p->m_id = particles().size();

    // Particles without production vertex are added to the root vertex
    if ( !p->production_vertex() ) {
        m_rootvertex->add_particle_out(p);
    }
}


GenEvent::GenEvent(const GenEvent&e) {
    if (this != &e)
    {
        std::lock(m_lock_attributes, e.m_lock_attributes);
        std::lock_guard<std::recursive_mutex> lhs_lk(m_lock_attributes, std::adopt_lock);
        std::lock_guard<std::recursive_mutex> rhs_lk(e.m_lock_attributes, std::adopt_lock);
        GenEventData tdata;
        e.write_data(tdata);
        read_data(tdata);
        m_run_info = e.m_run_info;
    }
}

GenEvent::~GenEvent() {
    for ( auto attm = m_attributes.begin(); attm != m_attributes.end(); ++attm) {
        for ( auto att = attm->second.begin(); att != attm->second.end(); ++att) { if (att->second) att->second->m_event = nullptr;}
    }
    for  ( auto v = m_vertices.begin(); v != m_vertices.end(); ++v ) if (*v)  if ((*v)->m_event == this) (*v)->m_event = nullptr;
    for  ( auto p = m_particles.begin(); p != m_particles.end(); ++p ) if (*p)  if ((*p)->m_event == this)  (*p)->m_event = nullptr;
}

GenEvent& GenEvent::operator=(const GenEvent& e) {
    if (this != &e)
    {
        std::lock(m_lock_attributes, e.m_lock_attributes);
        std::lock_guard<std::recursive_mutex> lhs_lk(m_lock_attributes, std::adopt_lock);
        std::lock_guard<std::recursive_mutex> rhs_lk(e.m_lock_attributes, std::adopt_lock);
        GenEventData tdata;
        e.write_data(tdata);
        read_data(tdata);
        m_run_info = e.m_run_info;
    }
    return *this;
}


void GenEvent::add_vertex(GenVertexPtr v) {
    if ( !v|| v->in_event() ) return;
    m_vertices.emplace_back(v);

    v->m_event = this;
    v->m_id = -(int)vertices().size();

    // Add all incoming and outgoing particles and restore their production/end vertices
    for (const auto& p: v->particles_in()) {
        if (!p->in_event()) add_particle(p);
        p->m_end_vertex = v->shared_from_this();
    }

    for (const auto& p: v->particles_out()) {
        if (!p->in_event()) add_particle(p);
        p->m_production_vertex = v;
    }
}


void GenEvent::remove_particle(GenParticlePtr p) {
    if ( !p || p->parent_event() != this ) return;

    HEPMC3_DEBUG(30, "GenEvent::remove_particle - called with particle: " << p->id());
    GenVertexPtr end_vtx = p->end_vertex();
    if ( end_vtx ) {
        end_vtx->remove_particle_in(p);

        // If that was the only incoming particle, remove vertex from the event
        if ( end_vtx->particles_in().empty() )  remove_vertex(end_vtx);
    }

    GenVertexPtr prod_vtx = p->production_vertex();
    if ( prod_vtx ) {
        prod_vtx->remove_particle_out(p);

        // If that was the only outgoing particle, remove vertex from the event
        if ( prod_vtx->particles_out().empty() ) remove_vertex(prod_vtx);
    }

    HEPMC3_DEBUG(30, "GenEvent::remove_particle - erasing particle: " << p->id())

    int idx = p->id();
    auto it = m_particles.erase(m_particles.begin() + idx-1);

    // Remove attributes of this particle
    std::lock_guard<std::recursive_mutex> lock(m_lock_attributes);
    for (att_key_t& vt1: m_attributes) {
        auto vt2 = vt1.second.find(idx);
        if (vt2 == vt1.second.end()) continue;
        vt1.second.erase(vt2);
    }

    //
    // Reassign id of attributes with id above this one
    //
    std::vector< std::pair< int, std::shared_ptr<Attribute> > > changed_attributes;

    for (att_key_t& vt1: m_attributes) {
        changed_attributes.clear();

        for (auto vt2 = vt1.second.begin(); vt2 != vt1.second.end(); ++vt2) {
            if ( (*vt2).first > p->id() ) {
                changed_attributes.emplace_back(*vt2);
            }
        }

        std::sort(changed_attributes.begin(),changed_attributes.end(), [](const std::pair< int, std::shared_ptr<Attribute> > &a, const std::pair< int, std::shared_ptr<Attribute> > &b) { return a.first < b.first; });
        for ( const auto& val: changed_attributes ) {
            vt1.second.erase(val.first);
            vt1.second[val.first-1] = val.second;
        }
    }
    // Reassign id of particles with id above this one
    for (; it != m_particles.end(); ++it) {
        --((*it)->m_id);
    }

    // Finally - set parent event and id of this particle to 0
    p->m_event = nullptr;
    p->m_id    = 0;
}

void GenEvent::remove_particles(std::vector<GenParticlePtr> v) {
    std::sort(v.begin(), v.end(), [](const GenParticlePtr& p1, const GenParticlePtr& p2) { return p1->id() > p2->id();});

    for (auto p = v.begin(); p != v.end(); ++p) {
        remove_particle(*p);
    }
}

void GenEvent::remove_vertex(GenVertexPtr v) {
    if ( !v || v->parent_event() != this ) return;

    HEPMC3_DEBUG(30, "GenEvent::remove_vertex   - called with vertex:  " << v->id());
    std::shared_ptr<GenVertex> null_vtx;

    for (const auto& p: v->particles_in()) {
        p->m_end_vertex = std::weak_ptr<GenVertex>();
    }

    for (const auto& p: v->particles_out()) {
        p->m_production_vertex = std::weak_ptr<GenVertex>();

        // recursive delete rest of the tree
        remove_particle(p);
    }

    // Erase this vertex from vertices list
    HEPMC3_DEBUG(30, "GenEvent::remove_vertex   - erasing vertex: " << v->id())

    int idx = -v->id();
    auto it = m_vertices.erase(m_vertices.begin() + idx-1);
    // Remove attributes of this vertex
    std::lock_guard<std::recursive_mutex> lock(m_lock_attributes);
    for (att_key_t& vt1: m_attributes) {
        auto vt2 = vt1.second.find(-idx);
        if (vt2 == vt1.second.end()) continue;
        vt1.second.erase(vt2);
    }

    //
    // Reassign id of attributes with id below this one
    //

    std::vector< std::pair< int, std::shared_ptr<Attribute> > > changed_attributes;

    for ( att_key_t& vt1: m_attributes ) {
        changed_attributes.clear();

        for (auto vt2 = vt1.second.begin(); vt2 != vt1.second.end(); ++vt2) {
            if ( (*vt2).first < v->id() ) {
                changed_attributes.emplace_back(*vt2);
            }
        }

        std::reverse(changed_attributes.begin(),changed_attributes.end());
        std::sort(changed_attributes.begin(),changed_attributes.end(),[](const std::pair< int, std::shared_ptr<Attribute> > &a, const std::pair< int, std::shared_ptr<Attribute> > &b) { return a.first > b.first; });
        for ( const auto& val: changed_attributes ) {
            vt1.second.erase(val.first);
            vt1.second[val.first+1] = val.second;
        }
    }

    // Reassign id of particles with id above this one
    for (; it != m_vertices.end(); ++it) {
        ++((*it)->m_id);
    }

    // Finally - set parent event and id of this vertex to 0
    v->m_event = nullptr;
    v->m_id    = 0;
}
/* This looks dangerously similar to the recusive event traversel that we forbade in the
       Core library due to wories about generator dependence
*/
static bool visit_children(std::map<ConstGenVertexPtr, int>  &a, const ConstGenVertexPtr& v)
{
    for (const ConstGenParticlePtr& p: v->particles_out()) {
        if (p->end_vertex())
        {
            if (a[p->end_vertex()] != 0) { return true; }
            a[p->end_vertex()]++;
            if (visit_children(a, p->end_vertex())) return true;
        }
    }
    return false;
}

void GenEvent::add_tree(const std::vector<GenParticlePtr> &parts) {
    m_particles.reserve(m_particles.size() + parts.size());
    m_vertices.reserve(m_vertices.size() + parts.size());
    std::shared_ptr<IntAttribute> existing_hc = attribute<IntAttribute>("cycles");
    bool has_cycles = false;
    std::map<GenVertexPtr, int>  sortingv;
    std::vector<GenVertexPtr> noinv;
    if (existing_hc)     if (existing_hc->value() != 0) has_cycles = true;
    if (!existing_hc)
    {
        for (const GenParticlePtr& p: parts) {
            GenVertexPtr v = p->production_vertex();
            if (v) sortingv[v]=0;
            if ( !v || v->particles_in().empty()) {
                GenVertexPtr v2 = p->end_vertex();
                if (v2) {noinv.emplace_back(v2); sortingv[v2] = 0;}
            }
        }
        for (const GenVertexPtr& v: noinv) {
            std::map<ConstGenVertexPtr, int>  sorting_temp(sortingv.begin(), sortingv.end());
            has_cycles = (has_cycles || visit_children(sorting_temp, v));
        }
    }
    if (has_cycles) {
        add_attribute("cycles", std::make_shared<IntAttribute>(1));
        /* Commented out  as improvemnts allow us to do sorting in other way.
         for ( std::map<GenVertexPtr,int>::iterator vi=sortingv.begin();vi!=sortingv.end();++vi) if ( !vi->first->in_event() ) add_vertex(vi->first);
         return;
         */
    }

    std::deque<GenVertexPtr> sorting;

    // Find all starting vertices (end vertex of particles that have no production vertex)
    for (const auto& p: parts) {
        const GenVertexPtr &v = p->production_vertex();
        if ( !v || v->particles_in().empty() ) {
            const GenVertexPtr &v2 = p->end_vertex();
            if (v2) sorting.emplace_back(v2);
        }
    }

    HEPMC3_DEBUG_CODE_BLOCK(
        unsigned int sorting_loop_count = 0;
        unsigned int max_deque_size     = 0;
    )

    // Add vertices to the event in topological order
    while ( !sorting.empty() ) {
        HEPMC3_DEBUG_CODE_BLOCK(
            if ( sorting.size() > max_deque_size ) max_deque_size = sorting.size();
            ++sorting_loop_count;
        )

            GenVertexPtr &v = sorting.front();

        bool added = false;

        // Add all mothers to the front of the list
        for (const auto& p: v->particles_in() ) {
            GenVertexPtr v2 = p->production_vertex();
            if ( v2 && !v2->in_event() && find(sorting.begin(), sorting.end(), v2) == sorting.end() ) {
                sorting.push_front(v2);
                added = true;
            }
        }

        // If we have added at least one production vertex,
        // our vertex is not the first one on the list
        if ( added ) continue;

        // If vertex not yet added
        if ( !v->in_event() ) {
            add_vertex(v);

            // Add all end vertices to the end of the list
            for (const auto& p: v->particles_out()) {
                GenVertexPtr v2 = p->end_vertex();
                if ( v2 && !v2->in_event()&& find(sorting.begin(), sorting.end(), v2) == sorting.end() ) {
                    sorting.emplace_back(v2);
                }
            }
        }

        sorting.pop_front();
    }

    // LL: Make sure root vertex has index zero and is not written out
    if ( m_rootvertex->id() != 0 ) {
        const int vx = -1 - m_rootvertex->id();
        const int rootid = m_rootvertex->id();
        if ( vx >= 0 && vx < (int) m_vertices.size() && m_vertices[vx] == m_rootvertex ) {
            auto next = m_vertices.erase(m_vertices.begin() + vx);
            std::lock_guard<std::recursive_mutex> lock(m_lock_attributes);
            for (auto & vt1: m_attributes) {
                std::vector< std::pair< int, std::shared_ptr<Attribute> > > changed_attributes;
                for ( const auto& vt2 : vt1.second ) {
                    if ( vt2.first <= rootid ) {
                        changed_attributes.emplace_back(vt2);
                    }
                }
                for ( const auto& val : changed_attributes ) {
                    vt1.second.erase(val.first);
                    vt1.second[val.first == rootid? 0: val.first + 1] = val.second;
                }
            }
            m_rootvertex->m_id = 0;
            while ( next != m_vertices.end() ) {
                ++((*next++)->m_id);
            }
        } else {
            HEPMC3_WARNING_LEVEL(700,"GenEvent::add_tree Suspicious looking rootvertex found. Will try to cope.")
        }
    }

    HEPMC3_DEBUG_CODE_BLOCK(
        HEPMC3_DEBUG(6, "GenEvent - particles sorted: "
                     << this->particles().size() << ", max deque size: "
                     << max_deque_size << ", iterations: " << sorting_loop_count)
    )
}


void GenEvent::reserve(const size_t& parts, const size_t& verts) {
    m_particles.reserve(parts);
    m_vertices.reserve(verts);
}


void GenEvent::set_units(Units::MomentumUnit new_momentum_unit, Units::LengthUnit new_length_unit) {
    if ( new_momentum_unit != m_momentum_unit ) {
        for ( GenParticlePtr& p: m_particles ) {
            Units::convert(p->m_data.momentum, m_momentum_unit, new_momentum_unit);
            Units::convert(p->m_data.mass, m_momentum_unit, new_momentum_unit);
        }

        m_momentum_unit = new_momentum_unit;
    }

    if ( new_length_unit != m_length_unit ) {
        for (GenVertexPtr& v: m_vertices) {
            FourVector &fv = v->m_data.position;
            if ( !fv.is_zero() ) Units::convert( fv, m_length_unit, new_length_unit );
        }

        m_length_unit = new_length_unit;
    }
}


const FourVector& GenEvent::event_pos() const {
    return m_rootvertex->data().position;
}

std::vector<ConstGenParticlePtr> GenEvent::beams(const int status) const {
    if (!status) return std::const_pointer_cast<const GenVertex>(m_rootvertex)->particles_out();
    std::vector<ConstGenParticlePtr> ret;
    for (auto& p: m_rootvertex->particles_out()) if (p->status() == status) ret.emplace_back(p);
    return ret;
}

std::vector<ConstGenParticlePtr> GenEvent::beams() const {
    return std::const_pointer_cast<const GenVertex>(m_rootvertex)->particles_out();
}


const std::vector<GenParticlePtr> & GenEvent::beams() {
    return m_rootvertex->particles_out();
}

void GenEvent::shift_position_by(const FourVector & delta) {
    m_rootvertex->set_position(event_pos() + delta);

    // Offset all vertices
    for ( GenVertexPtr& v: m_vertices ) {
        if ( v->has_set_position() ) {
            v->set_position(v->position() + delta);
        }
    }
}

bool GenEvent::rotate(const FourVector&  delta)
{
    long double cosa = std::cos(delta.x());
    long double sina = std::sin(delta.x());
    long double cosb = std::cos(delta.y());
    long double sinb = std::sin(delta.y());
    long double cosg = std::cos(delta.z());
    long double sing = std::sin(delta.z());

    for ( auto& p: m_particles)
    {
        const FourVector& mom = p->momentum();
        long double tempX = mom.x();
        long double tempY = mom.y();
        long double tempZ = mom.z();

        long double tempY_ = cosa*tempY+sina*tempZ;
        long double tempZ_ = -sina*tempY+cosa*tempZ;
        tempY = tempY_;
        tempZ = tempZ_;

        long double tempX_ = cosb*tempX-sinb*tempZ;
        tempZ_ = sinb*tempX+cosb*tempZ;
        tempX = tempX_;
        tempZ = tempZ_;

        tempX_ = cosg*tempX+sing*tempY;
        tempY_ = -sing*tempX+cosg*tempY;
        tempX = tempX_;
        tempY = tempY_;

        FourVector temp(tempX, tempY, tempZ, mom.e());
        p->set_momentum(temp);
    }
    for (auto& v: m_vertices)
    {
        const FourVector& pos = v->position();
        if (pos.is_zero()) continue;

        long double tempX = pos.x();
        long double tempY = pos.y();
        long double tempZ = pos.z();

        long double tempY_ = cosa*tempY+sina*tempZ;
        long double tempZ_ = -sina*tempY+cosa*tempZ;
        tempY = tempY_;
        tempZ = tempZ_;

        long double tempX_ = cosb*tempX-sinb*tempZ;
        tempZ_ = sinb*tempX+cosb*tempZ;
        tempX = tempX_;
        tempZ = tempZ_;

        tempX_ = cosg*tempX+sing*tempY;
        tempY_ = -sing*tempX+cosg*tempY;
        tempX = tempX_;
        tempY = tempY_;

        FourVector temp(tempX, tempY, tempZ, pos.t());
        v->set_position(temp);
    }


    return true;
}

bool GenEvent::reflect(const int axis)
{
    if ( axis > 3 || axis < 0 )
    {
        HEPMC3_WARNING_LEVEL(400,"GenEvent::reflect: wrong axis")
        return false;
    }
    switch (axis)
    {
    case 0:
        for ( auto& p: m_particles) { FourVector temp = p->momentum(); temp.setX(-p->momentum().x()); p->set_momentum(temp);}
        for ( auto& v: m_vertices)  { FourVector temp = v->position(); temp.setX(-v->position().x()); v->set_position(temp);}
        break;
    case 1:
        for ( auto& p: m_particles) { FourVector temp = p->momentum(); temp.setY(-p->momentum().y()); p->set_momentum(temp);}
        for ( auto& v: m_vertices)  { FourVector temp = v->position(); temp.setY(-v->position().y()); v->set_position(temp);}
        break;
    case 2:
        for ( auto& p: m_particles) { FourVector temp = p->momentum(); temp.setZ(-p->momentum().z()); p->set_momentum(temp);}
        for ( auto& v: m_vertices)  { FourVector temp = v->position(); temp.setZ(-v->position().z()); v->set_position(temp);}
        break;
    case 3:
        for ( auto& p: m_particles) { FourVector temp = p->momentum(); temp.setT(-p->momentum().e()); p->set_momentum(temp);}
        for ( auto& v: m_vertices)  { FourVector temp = v->position(); temp.setT(-v->position().t()); v->set_position(temp);}
        break;
    default:
        return false;
    }

    return true;
}

bool GenEvent::boost(const FourVector&  delta)
{
    double deltalength2 = delta.length2();
    if (deltalength2 > 1.0)
    {
        HEPMC3_WARNING_LEVEL(400,"GenEvent::boost: wrong large boost vector. Will leave event as is.")
        return false;
    }
    if (std::abs(deltalength2-1.0) < std::numeric_limits<double>::epsilon())
    {
        HEPMC3_WARNING_LEVEL(400,"GenEvent::boost: too large gamma. Will leave event as is.")
        return false;
    }
    if (std::abs(deltalength2) < std::numeric_limits<double>::epsilon())
    {
        HEPMC3_WARNING_LEVEL(400,"GenEvent::boost: wrong small boost vector. Will leave event as is.")
        return true;
    }
    long double deltaX = delta.x();
    long double deltaY = delta.y();
    long double deltaZ = delta.z();
    long double deltalength = std::sqrt(deltalength2);
    long double gamma = 1.0/std::sqrt(1.0-deltalength2);

    for ( auto& p: m_particles)
    {
        const FourVector& mom = p->momentum();

        long double tempX = mom.x();
        long double tempY = mom.y();
        long double tempZ = mom.z();
        long double tempE = mom.e();
        long double nr = (deltaX*tempX+deltaY*tempY+deltaZ*tempZ)/deltalength;
        long double gfac = (gamma-1)*nr/deltalength-tempE*gamma;
        tempX+=(deltaX*gfac);
        tempY+=(deltaY*gfac);
        tempZ+=(deltaZ*gfac);
        tempE = gamma*(tempE-deltalength*nr);
        FourVector temp(tempX, tempY, tempZ, tempE);
        p->set_momentum(temp);
    }

    return true;
}

void GenEvent::clear() {
    std::lock_guard<std::recursive_mutex> lock(m_lock_attributes);
    m_event_number = 0;
    m_rootvertex = std::make_shared<GenVertex>();
    m_weights.clear();
    m_attributes.clear();
    m_particles.clear();
    m_vertices.clear();
}

void GenEvent::remove_attribute(const std::string &name,  const int& id) {
    std::lock_guard<std::recursive_mutex> lock(m_lock_attributes);
    auto i1 = m_attributes.find(name);
    if ( i1 == m_attributes.end() ) return;

    auto i2 = i1->second.find(id);
    if ( i2 == i1->second.end() ) return;

    i1->second.erase(i2);
}

std::vector<std::string> GenEvent::attribute_names(const int& id) const {
    std::vector<std::string> results;

    for (const att_key_t& vt1: m_attributes) {
        if ( vt1.second.count(id) == 1 ) {
            results.emplace_back(vt1.first);
        }
    }

    return results;
}

void GenEvent::write_data(GenEventData& data) const {
    // Reserve memory for containers
    data.particles.reserve(this->particles().size());
    data.vertices.reserve(this->vertices().size());
    data.links1.reserve(this->particles().size()*2);
    data.links2.reserve(this->particles().size()*2);
    data.attribute_id.reserve(m_attributes.size());
    data.attribute_name.reserve(m_attributes.size());
    data.attribute_string.reserve(m_attributes.size());

    // Fill event data
    data.event_number  = this->event_number();
    data.momentum_unit = this->momentum_unit();
    data.length_unit   = this->length_unit();
    data.event_pos     = this->event_pos();

    // Fill containers
    data.weights = this->weights();

    for (const ConstGenParticlePtr& p: this->particles()) {
        data.particles.emplace_back(p->data());
    }

    for (const ConstGenVertexPtr& v: this->vertices()) {
        data.vertices.emplace_back(v->data());
        int v_id = v->id();

        for (const ConstGenParticlePtr& p: v->particles_in()) {
            data.links1.emplace_back(p->id());
            data.links2.emplace_back(v_id);
        }

        for (const ConstGenParticlePtr& p: v->particles_out()) {
            data.links1.emplace_back(v_id);
            data.links2.emplace_back(p->id());
        }
    }

    for (const att_key_t& vt1: this->attributes()) {
        for (const att_val_t& vt2: vt1.second) {
            std::string st;

            bool status = vt2.second->to_string(st);

            if ( !status ) {
                HEPMC3_WARNING_LEVEL(300,"GenEvent::write_data: problem serializing attribute: " << vt1.first)
            }
            else {
                data.attribute_id.emplace_back(vt2.first);
                data.attribute_name.emplace_back(vt1.first);
                data.attribute_string.emplace_back(st);
            }
        }
    }
}


void GenEvent::read_data(const GenEventData &data) {
    this->clear();
    this->set_event_number(data.event_number);
    //Note: set_units checks the current unit of event, i.e. applicable only for fully constructed event.
    m_momentum_unit = data.momentum_unit;
    m_length_unit = data.length_unit;
    this->shift_position_to(data.event_pos);

    // Fill weights
    this->weights() = data.weights;
    m_particles.reserve(data.particles.size());
    m_vertices.reserve(data.vertices.size());

    // Fill particle information
    for ( const GenParticleData &pd: data.particles ) {
        m_particles.emplace_back(std::make_shared<GenParticle>(pd));
        m_particles.back()->m_event = this;
        m_particles.back()->m_id    = m_particles.size();
    }

    // Fill vertex information
    for ( const GenVertexData &vd: data.vertices ) {
        m_vertices.emplace_back(std::make_shared<GenVertex>(vd));
        m_vertices.back()->m_event = this;
        m_vertices.back()->m_id    = -(int)m_vertices.size();
    }

    // Restore links
    for (unsigned int i = 0; i < data.links1.size(); ++i) {
        const int id1 = data.links1[i];
        const int id2 = data.links2[i];
        /* @note:
        The  meaningfull combinations for (id1,id2) are:
        (+-)  --  particle has end vertex
        (-+)  --  particle  has production vertex
        */
        if ((id1 < 0 && id2 <0) || (id1 > 0 && id2 > 0))   {
            HEPMC3_WARNING_LEVEL(600,"GenEvent::read_data: wrong link: " << id1 << " " << id2);
            continue;
        }

        if ( id1 > 0 ) { m_vertices[ (-id2)-1 ]->add_particle_in ( m_particles[ id1-1 ] ); continue; }
        if ( id1 < 0 ) { m_vertices[ (-id1)-1 ]->add_particle_out( m_particles[ id2-1 ] );   continue; }
    }
    for (auto& p:  m_particles) if (!p->production_vertex()) m_rootvertex->add_particle_out(p);

    // Read attributes
    std::lock_guard<std::recursive_mutex> lock(m_lock_attributes);
    for (unsigned int i = 0; i < data.attribute_id.size(); ++i) {
        ///Disallow empty strings
        const std::string name = data.attribute_name[i];
        if (name.length() == 0) continue;
        const int id = data.attribute_id[i];
        if (m_attributes.count(name) == 0) m_attributes[name] = std::map<int, std::shared_ptr<Attribute> >();
        auto att = std::make_shared<StringAttribute>(data.attribute_string[i]);
        att->m_event = this;
        if ( id > 0 && id <= int(m_particles.size()) ) {
            att->m_particle = m_particles[id - 1];
        }
        if ( id < 0 && -id <= int(m_vertices.size()) ) {
            att->m_vertex = m_vertices[-id - 1];
        }
        m_attributes[name][id] = att;
    }
}


//
// Deprecated functions
//

void GenEvent::set_beam_particles(GenParticlePtr p1, GenParticlePtr p2) {
    m_rootvertex->add_particle_out(p1);
    m_rootvertex->add_particle_out(p2);
}

void GenEvent::add_beam_particle(GenParticlePtr p1) {
    if (!p1)
    {
        HEPMC3_WARNING_LEVEL(700,"Attempting to add an empty particle as beam particle. Ignored.")
        return;
    }
    if (p1->in_event() && p1->parent_event() != this)
    {
        HEPMC3_WARNING_LEVEL(700,"Attempting to add particle from another event. Ignored.")
        return;
    }
    if (p1->production_vertex())  p1->production_vertex()->remove_particle_out(p1);
    //Particle w/o production vertex is added to root vertex.
    add_particle(p1);
    p1->set_status(4);
}


std::string GenEvent::attribute_as_string(const std::string &name, const int& id) const {
    std::lock_guard<std::recursive_mutex> lock(m_lock_attributes);
    auto i1 = m_attributes.find(name);
    if ( i1 == m_attributes.end() ) {
        if ( id == 0 && run_info() ) {
            return run_info()->attribute_as_string(name);
        }
        return {};
    }

    auto i2 = i1->second.find(id);
    if (i2 == i1->second.end() ) return {};

    if ( !i2->second ) return {};

    std::string ret;
    i2->second->to_string(ret);

    return ret;
}

void GenEvent::add_attribute(const std::string &name, const std::shared_ptr<Attribute> &att, const int& id ) {
    ///Disallow empty strings
    if (name.length() == 0) return;
    if (!att)  return;
    std::lock_guard<std::recursive_mutex> lock(m_lock_attributes);
    if (m_attributes.count(name) == 0) m_attributes[name] = std::map<int, std::shared_ptr<Attribute> >();
    m_attributes[name][id] = att;
    att->m_event = this;
    if ( id > 0 && id <= int(particles().size()) ) {
        att->m_particle = particles()[id - 1];
    }
    if ( id < 0 && -id <= int(vertices().size()) ) {
        att->m_vertex = vertices()[-id - 1];
    }
}


void GenEvent::add_attributes(const std::vector<std::string> &names, const std::vector<std::shared_ptr<Attribute> > &atts, const std::vector<int>& ids) {
    size_t N = names.size();
    if ( N == 0 ) return;
    if (N != atts.size()) return;
    if (N != ids.size()) return;

    std::vector<std::string> unames = names;
    vector<std::string>::iterator ip;
    ip = std::unique(unames.begin(), unames.end());
    unames.resize(std::distance(unames.begin(), ip));
    std::lock_guard<std::recursive_mutex> lock(m_lock_attributes);
    for (const auto& name: unames) {
        if (m_attributes.count(name) == 0) m_attributes[name] = std::map<int, std::shared_ptr<Attribute> >();
    }
    const int particles_size = int(m_particles.size());
    const int vertices_size = int(m_vertices.size());
    for (size_t i = 0; i < N; i++) {
        ///Disallow empty strings
        if (names.at(i).length() == 0) continue;
        if (!atts[i])  continue;
        m_attributes[names.at(i)][ids.at(i)] = atts[i];
        atts[i]->m_event = this;
        if ( ids.at(i) > 0 && ids.at(i) <= particles_size )
        { atts[i]->m_particle = m_particles[ids.at(i) - 1]; }
        else {
            if ( ids.at(i) < 0 && -ids.at(i) <= vertices_size ) {
                atts[i]->m_vertex = m_vertices[-ids.at(i) - 1];
            }
        }
    }
}

void GenEvent::add_attributes(const std::string& name, const std::vector<std::shared_ptr<Attribute> > &atts, const std::vector<int>& ids) {
    if (name.length() == 0) return;
    size_t N = ids.size();
    if(!N) return;
    if ( N != atts.size()) return;

    std::lock_guard<std::recursive_mutex> lock(m_lock_attributes);
    if (m_attributes.count(name) == 0) m_attributes[name] = std::map<int, std::shared_ptr<Attribute> >();
    auto& tmap = m_attributes[name];
    const int particles_size = int(m_particles.size());
    const int vertices_size = int(m_vertices.size());
    for (size_t i = 0; i < N; i++) {
        ///Disallow empty strings
        if (!atts[i])  continue;
        tmap[ids.at(i)] = atts[i];
        atts[i]->m_event = this;
        if ( ids.at(i) > 0 && ids.at(i) <= particles_size )
        { atts[i]->m_particle = m_particles[ids.at(i) - 1]; }
        else {
            if ( ids.at(i) < 0 && -ids.at(i) <= vertices_size ) {
                atts[i]->m_vertex = m_vertices[-ids.at(i) - 1];
            }
        }
    }
}
void GenEvent::add_attributes(const std::string& name, const std::vector<std::pair<int, std::shared_ptr<Attribute> > > &atts) {
    if (name.length() == 0) return;
    if (atts.empty()) return;
    std::lock_guard<std::recursive_mutex> lock(m_lock_attributes);
    if (m_attributes.count(name) == 0) m_attributes[name] = std::map<int, std::shared_ptr<Attribute> >();
    auto& tmap = m_attributes[name];
    const int particles_size = int(m_particles.size());
    const int vertices_size = int(m_vertices.size());
    for (const auto& att: atts) {
        ///Disallow empty strings
        if (!att.second)  continue;
        tmap.insert(att);
        att.second->m_event = this;
        if ( att.first > 0 && att.first <= particles_size )
        { att.second->m_particle = m_particles[att.first - 1]; }
        else {
            if ( att.first < 0 && -att.first <= vertices_size ) {
                att.second->m_vertex = m_vertices[-att.first - 1];
            }
        }
    }
}

} // namespace HepMC3

// ============================================================
// hepmc3/Print.cc
// Copyright (C) 2014-2023 The HepMC collaboration
// ============================================================
// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2023 The HepMC collaboration (see AUTHORS for details)
//
///
/// @file Print.cc
/// @brief Implementation of static \b class Print
///
///
#include "HepMC3/Print.h"
#include "HepMC3/Attribute.h"


namespace HepMC3 {

void Print::content(std::ostream& os, const GenEvent &event) {
    os << "--------------------------------" << std::endl;
    os << "--------- EVENT CONTENT --------" << std::endl;
    os << "--------------------------------" << std::endl;
    os << std::endl;

    os << "Weights (" << event.weights().size() << "): " << std::endl;
    for (const auto& w: event.weights()) {
        os << " " << w;
    }

    os << "Attributes:" << std::endl;

    for (const auto& vt1: event.attributes()) {
        for (const auto& vt2: vt1.second) {
            os << vt2.first << ": " << vt1.first << std::endl;
        }
    }

    os << "GenParticlePtr (" << event.particles().size() << ")" << std::endl;

    for (const ConstGenParticlePtr& p: event.particles()) {
        Print::line(os, p, true);
        os << std::endl;
    }

    os << "GenVertexPtr (" << event.vertices().size() << ")" << std::endl;
    for ( const ConstGenVertexPtr& v: event.vertices() ) {
        Print::line(os, v, true);
        os << std::endl;
    }

    os << "-----------------------------" << std::endl;
}

void Print::listing(std::ostream& os, const GenEvent &event, unsigned short precision) {
    // Find the current stream state
    std::ios_base::fmtflags orig = os.flags();
    std::streamsize         prec = os.precision();

    // Set precision
    os.precision(precision);

    os << "________________________________________________________________________" << std::endl;
    os << "GenEvent: #" << event.event_number() << std::endl;
    os << " Momentum units: " << Units::name(event.momentum_unit())
       << " Position units: " << Units::name(event.length_unit()) << std::endl;
    os << " Entries in this event: " << event.vertices().size() << " vertices, "
       << event.particles().size() << " particles, "
       << event.weights().size()   << " weights." << std::endl;

    const FourVector &pos = event.event_pos();
    os << " Position offset: " << pos.x() << ", " << pos.y() << ", " << pos.z() << ", " << pos.t() << std::endl;

    // Print a legend to describe the particle info
    os << "                                    GenParticle Legend" << std::endl;
    os << "         ID    PDG ID   "
       << "( px,       py,       pz,     E )"
       << "   Stat ProdVtx" << std::endl;
    os << "________________________________________________________________________" << std::endl;

    // Print all vertices
    for (const ConstGenVertexPtr& v: event.vertices()) {
        Print::listing(os, v);
    }

    // Restore the stream state
    os.flags(orig);
    os.precision(prec);
    os << "________________________________________________________________________" << std::endl;
}

void Print::listing(std::ostream& os, const GenRunInfo &ri, unsigned short precision) {
    // Find the current stream state
    std::ios_base::fmtflags orig = os.flags();
    std::streamsize         prec = os.precision();

    // Set precision
    os.precision(precision);

    os << "________________________________________________________________________" << std::endl;
    os << "GenRunInfo:" << std::endl;

    std::vector<std::string> names = ri.weight_names();
    os << " Names: ( ";
    for (const auto& n: names) os << n;
    os << " )" << std::endl;

    os << " Tools: " << std::endl;

    for (const auto& t: ri.tools()) {
        Print::line(os, t);
    }
    os << "Attributes:" << std::endl;
    for (const auto& att: ri.attributes()) {
        std::string st;
        if ( !att.second->to_string(st) ) {
            HEPMC3_WARNING_LEVEL(300,"Print::listing: problem serializing attribute: " << att.first)
        }
        else { os << att.first << " " << st;}
        os << std::endl;
    }

    // Restore the stream state
    os.flags(orig);
    os.precision(prec);
    os << "________________________________________________________________________" << std::endl;
}

void Print::listing(std::ostream& os, ConstGenVertexPtr v) {
    if (!v) { os << "Vtx: Empty vertex" << std::endl; return;}
    os << "Vtx: ";
    os.width(6);
    os << v->id() << " stat: ";
    os.width(3);
    os << v->status();

    const FourVector &pos = v->position();
    if ( !pos.is_zero() ) {
        os << " (X,cT): " << pos.x() << " " << pos.y() << " " << pos.z() << " " << pos.t();
    }
    else os << " (X,cT): 0";

    os << std::endl;

    bool printed_header = false;

    // Print out all the incoming particles
    for (const ConstGenParticlePtr& p: v->particles_in()) {
        if ( !printed_header ) {
            os << " I: ";
            printed_header = true;
        }
        else os << "    ";

        Print::listing(os, p);
    }

    printed_header = false;

    // Print out all the outgoing particles
    for (const ConstGenParticlePtr& p: v->particles_out()) {
        if ( !printed_header ) {
            os << " O: ";
            printed_header = true;
        }
        else os << "    ";

        Print::listing(os, p);
    }
}

void Print::listing(std::ostream& os, ConstGenParticlePtr p) {
    if (!p) { os << " Empty particle" << std::endl; return;}
    os << " ";
    os.width(6);
    os << p->id();
    os.width(9);
    os << p->pid() << " ";
    os.width(9);
    os.setf(std::ios::scientific, std::ios::floatfield);
    os.setf(std::ios_base::showpos);

    const FourVector &momentum = p->momentum();

    os.width(9);
    os << momentum.px() << ",";
    os.width(9);
    os << momentum.py() << ",";
    os.width(9);
    os << momentum.pz() << ",";
    os.width(9);
    os << momentum.e() << " ";
    os.setf(std::ios::fmtflags(0), std::ios::floatfield);
    os.unsetf(std::ios_base::showpos);
    os.width(3);
    os << p->status();

    ConstGenVertexPtr prod = p->production_vertex();

    if ( prod ) {
        os.width(6);
        os << prod->id();
    }

    os << std::endl;
}
void Print::line(std::ostream& os, const GenEvent &event, bool attributes) {
    os << "GenEvent: #" << event.event_number();
    if (attributes) {
        for (const std::string& s: event.attribute_names()) {
            os << " " << s << "=" <<event.attribute_as_string(s);
        }
    }
}

void Print::line(std::ostream& os, const GenRunInfo &RunInfo, bool attributes) {
    os <<"GenRunInfo: Number of tools:" << RunInfo.tools().size();

    if (attributes) {
        for (const std::string& s: RunInfo.attribute_names()) {
            os << " " << s << "=" << RunInfo.attribute_as_string(s);
        }
    }
}

void Print::line(std::ostream& os, const GenRunInfo::ToolInfo& t) {
    os << "GenRunInfo::ToolInfo " << t.name<< " " << t.version << " " << t.description;
}

template <class T>
void line_v(std::ostream& os, T v, bool attributes) {
    if (!v) { os << "GenVertex: Empty" << std::endl; return;}
    os << "GenVertex:  " << v->id() << " stat: ";
    os.width(3);
    os << v->status();
    os << " in: "  << v->particles_in().size();
    os.width(3);
    os << " out: " << v->particles_out().size();

    const FourVector &pos = v->position();
    os << " has_set_position: ";
    if ( v->has_set_position() ) { os << "true"; }
    else  { os << "false"; }

    os << " (X,cT): " << pos.x() << ", " <<pos.y() << ", " << pos.z() << ", " << pos.t();
    if (attributes)
    {
        auto names     = v->attribute_names();
        for (const auto& ss: names) {
            os << " " << ss << "=" << (*v).attribute_as_string(ss);
        }
    }
}
void Print::line(std::ostream& os, ConstGenVertexPtr v, bool attributes) { line_v(os,v,attributes); }
void Print::line(std::ostream& os, GenVertexPtr v, bool attributes) { line_v(os,v,attributes); }



void Print::line(std::ostream& os, const FourVector& p) {
    os << "FourVector: ";
    // Find the current stream state
    std::ios_base::fmtflags orig = os.flags();
    os.setf(std::ios::scientific, std::ios::floatfield);
    os.setf(std::ios_base::showpos);
    std::streamsize prec = os.precision();
    // Set precision
    os.precision(2);
    os << " (P,E)=" << p.x()
       << "," << p.y()
       << "," << p.z()
       << "," << p.e();

    // Restore the stream state
    os.flags(orig);
    os.precision(prec);
}

template <class T>
void line_p(std::ostream& os, T p, bool attributes) {
    if (!p) { os << "GenParticle: Empty" << std::endl; return;}
    os << "GenParticle: ";
    os.width(3);
    os << p->id() <<" PDGID: ";
    os.width(5);
    os << p->pid();

    // Find the current stream state
    std::ios_base::fmtflags orig = os.flags();

    os.setf(std::ios::scientific, std::ios::floatfield);
    os.setf(std::ios_base::showpos);
    std::streamsize prec = os.precision();

    // Set precision
    os.precision(2);

    const FourVector &momentum = p->momentum();

    os << " (P,E)=" << momentum.px()
       << "," << momentum.py()
       << "," << momentum.pz()
       << "," << momentum.e();

    // Restore the stream state
    os.flags(orig);
    os.precision(prec);

    const ConstGenVertexPtr prod = p->production_vertex();
    const ConstGenVertexPtr end  = p->end_vertex();
    int prod_vtx_id   = (prod) ? prod->id() : 0;
    int end_vtx_id    = (end)  ? end->id()  : 0;
    auto names        = p->attribute_names();

    os << " Stat: " << p->status()
       << " PV: " << prod_vtx_id
       << " EV: " << end_vtx_id
       << " Attr: " << names.size();

    if (attributes)
    {
        for (const auto& ss: names) {
            os << " " << ss << "=" << (*p).attribute_as_string(ss);
        }
    }
}

void Print::line(std::ostream& os, ConstGenParticlePtr p, bool attributes) { line_p(os,p,attributes); }
void Print::line(std::ostream& os, GenParticlePtr p, bool attributes) { line_p(os,p,attributes); }

void Print::line(std::ostream& os, std::shared_ptr<GenCrossSection> &cs) {
    if (!cs) {os << " GenCrossSection: Empty"; return;}
    os << " GenCrossSection: " << cs->xsec(0)
       << " " << cs->xsec_err(0)
       << " " << cs->get_accepted_events()
       << " " << cs->get_attempted_events();
}

void Print::line(std::ostream& os, std::shared_ptr<GenHeavyIon> &hi) {
    if (!hi) {os << " GenHeavyIon: Empty"; return;}
    os << " GenHeavyIon: " << hi->Ncoll_hard
       << " " << hi->Npart_proj
       << " " << hi->Npart_targ
       << " " << hi->Ncoll
       << " " << hi->spectator_neutrons
       << " " << hi->spectator_protons
       << " " << hi->N_Nwounded_collisions
       << " " << hi->Nwounded_N_collisions
       << " " << hi->Nwounded_Nwounded_collisions
       << " " << hi->impact_parameter
       << " " << hi->event_plane_angle
       << " " << hi->eccentricity
       << " " << hi->sigma_inel_NN;
}

void Print::line(std::ostream& os, std::shared_ptr<GenPdfInfo> &pi) {
    if (!pi) {os << " GenPdfInfo: Empty"; return;}
    os << " GenPdfInfo: " << pi->parton_id[0]
       << " " << pi->parton_id[1]
       << " " << pi->x[0]
       << " " << pi->x[1]
       << " " << pi->scale
       << " " << pi->xf[0]
       << " " << pi->xf[1]
       << " " << pi->pdf_id[0]
       << " " << pi->pdf_id[1];
}

} // namespace HepMC3

// ============================================================
// hepmc3/WriterAscii.cc
// Copyright (C) 2014-2023 The HepMC collaboration
// ============================================================
// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2023 The HepMC collaboration (see AUTHORS for details)
//
///
/// @file WriterAscii.cc
/// @brief Implementation of \b class WriterAscii
///

#include <algorithm>//min max for VS2017
#include <cstring>


#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Units.h"
#include "HepMC3/Version.h"
#include "HepMC3/WriterAscii.h"

namespace HepMC3 {


WriterAscii::WriterAscii(const std::string &filename, std::shared_ptr<GenRunInfo> run)
    : m_file(filename),
      m_stream(&m_file)
{
    set_run_info(run);
    if ( !m_file.is_open() ) {
        HEPMC3_ERROR_LEVEL(200,"WriterAscii: could not open output file: " << filename)
    } else {
        const std::string header = "HepMC::Version " + version() + "\nHepMC::Asciiv3-START_EVENT_LISTING\n";
        m_file.write(header.data(), header.length());
        if ( run_info() ) write_run_info();
    }
    m_float_printf_specifier = " %." + std::to_string(m_precision) + "e";
    m_particle_printf_specifier = "P %i %i %i"
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier + " %i\n";
    m_vertex_short_printf_specifier = "V %i %i [%s]\n";
    m_vertex_long_printf_specifier = "V %i %i [%s] @"+ m_float_printf_specifier + m_float_printf_specifier + m_float_printf_specifier + m_float_printf_specifier + "\n";
}


WriterAscii::WriterAscii(std::ostream &stream, std::shared_ptr<GenRunInfo> run)
    : m_stream(&stream)
{
    set_run_info(run);
    const std::string header = "HepMC::Version " + version() + "\nHepMC::Asciiv3-START_EVENT_LISTING\n";
    m_stream->write(header.data(), header.length());
    if ( run_info() ) write_run_info();
    m_float_printf_specifier = " %." + std::to_string(m_precision) + "e";
    m_particle_printf_specifier = "P %i %i %i"
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier + " %i\n";
    m_vertex_short_printf_specifier = "V %i %i [%s]\n";
    m_vertex_long_printf_specifier = "V %i %i [%s] @"+ m_float_printf_specifier + m_float_printf_specifier + m_float_printf_specifier + m_float_printf_specifier + "\n";
}

WriterAscii::WriterAscii(std::shared_ptr<std::ostream> s_stream, std::shared_ptr<GenRunInfo> run)
    : m_shared_stream(s_stream),
      m_stream(s_stream.get())
{
    set_run_info(run);
    const std::string header = "HepMC::Version " + version() + "\nHepMC::Asciiv3-START_EVENT_LISTING\n";
    m_stream->write(header.data(), header.length());
    if ( run_info() ) write_run_info();
    m_float_printf_specifier = " %." + std::to_string(m_precision) + "e";
    m_particle_printf_specifier = "P %i %i %i"
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier + " %i\n";
    m_vertex_short_printf_specifier = "V %i %i [%s]\n";
    m_vertex_long_printf_specifier = "V %i %i [%s] @"+ m_float_printf_specifier + m_float_printf_specifier + m_float_printf_specifier + m_float_printf_specifier + "\n";
}

WriterAscii::~WriterAscii() {
    close();
    delete[] m_buffer;
}


void WriterAscii::write_event(const GenEvent &evt) {
    allocate_buffer();
    if ( !m_buffer ) return;
    auto float_printf_specifier_option = m_options.find("float_printf_specifier");
    std::string  letter=(float_printf_specifier_option != m_options.end())?float_printf_specifier_option->second.substr(0,2):"e";
    if (letter != "e" && letter != "E" && letter != "G" && letter != "g" && letter != "f" && letter != "F" ) letter = "e";
    m_float_printf_specifier = " %." + std::to_string(m_precision) + letter;


    m_particle_printf_specifier = "P %i %i %i"
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier + " %i\n";
    m_vertex_short_printf_specifier = "V %i %i [%s]\n";
    m_vertex_long_printf_specifier = "V %i %i [%s] @"+ m_float_printf_specifier + m_float_printf_specifier + m_float_printf_specifier + m_float_printf_specifier + "\n";

    // Make sure nothing was left from previous event
    flush();

    if ( !run_info() ) {
        set_run_info(evt.run_info());
        write_run_info();
    } else {
        if ( evt.run_info() && (run_info() != evt.run_info()) ) {
            HEPMC3_WARNING_LEVEL(600,"WriterAscii::write_event: GenEvents contain different GenRunInfo objects from - only the first such object will be serialized.")
        }
    }

    // Write event info
    flush();
    std::string especifier =  "E " + std::to_string(evt.event_number()) + " "
                              + std::to_string(evt.vertices().size()) + " "
                              + std::to_string(evt.particles().size());
    // Write event position if not zero
    const FourVector &pos = evt.event_pos();
    if ( !pos.is_zero() ) {
        especifier += ( " @"  + m_float_printf_specifier + m_float_printf_specifier + m_float_printf_specifier + m_float_printf_specifier + "\n" );
        m_cursor += sprintf(m_cursor, especifier.c_str(), pos.x(), pos.y(), pos.z(), pos.t());
    } else {
        m_cursor += sprintf(m_cursor, "%s\n", especifier.c_str());
    }
    flush();

    // Write units
    m_cursor += sprintf(m_cursor, "U %s %s\n", Units::name(evt.momentum_unit()).c_str(), Units::name(evt.length_unit()).c_str());
    flush();

    // Write weight values if present
    if ( !evt.weights().empty() ) {
        m_cursor += sprintf(m_cursor, "W");
        for (const auto& w: evt.weights())
        {
            m_cursor += sprintf(m_cursor, " %.*e", std::min(3*m_precision, 22), w);
            flush();
        }
        m_cursor += sprintf(m_cursor, "\n");
        flush();
    }

    // Write attributes
    for ( const auto& vt1: evt.attributes() ) {
        for ( const auto& vt2: vt1.second ) {
            std::string st;
            bool status = vt2.second->to_string(st);

            if ( !status ) {
                HEPMC3_WARNING_LEVEL(300,"WriterAscii::write_event: problem serializing attribute: " << vt1.first)
            }
            else {
                m_cursor += sprintf(m_cursor, "A %i ", vt2.first);
                write_string(escape(vt1.first));
                flush();
                m_cursor += sprintf(m_cursor, " ");
                write_string(escape(st));
                m_cursor += sprintf(m_cursor, "\n");
                flush();
            }
        }
    }


    // Print particles
    std::map<int, bool> alreadywritten;
    for (const ConstGenParticlePtr& p: evt.particles()) {
        // Check to see if we need to write a vertex first
        ConstGenVertexPtr v = p->production_vertex();
        int parent_object = 0;

        if (v) {
            // Check if we need this vertex at all
            // Yes, use vertex as parent object
            if ( v->particles_in().size() > 1 || !v->data().is_zero() ) { parent_object = v->id(); }
            // No, use particle as parent object
            // Add check for attributes of this vertex
            else {
                if ( v->particles_in().size() == 1 )                  { parent_object = v->particles_in().front()->id();}
                else {if ( v->particles_in().empty() ) {HEPMC3_DEBUG(30, "WriterAscii::write_event - found a vertex without incoming particles: " << v->id());}}
            }
            // Usage of map instead of simple counter helps to deal with events with random ids of vertices.
            if (alreadywritten.count(v->id()) == 0 && parent_object < 0)
            { write_vertex(v); alreadywritten[v->id()] = true; }
        }

        write_particle(p, parent_object);
    }
    alreadywritten.clear();

    // Flush rest of the buffer to file
    forced_flush();
}


void WriterAscii::allocate_buffer() {
    if ( m_buffer ) return;
    while ( m_buffer == nullptr && m_buffer_size >= 512 ) {
        try {
            m_buffer = new char[ m_buffer_size ]();
        } catch (const std::bad_alloc& e) {
            delete[] m_buffer;
            m_buffer_size /= 2;
            HEPMC3_WARNING_LEVEL(200,"WriterAscii::allocate_buffer:" << e.what() << " buffer size too large. Dividing by 2. New size: " << m_buffer_size)
        }
    }

    if ( !m_buffer ) {
        HEPMC3_ERROR_LEVEL(200,"WriterAscii::allocate_buffer: could not allocate buffer!")
        return;
    }
    m_cursor = m_buffer;
}


std::string WriterAscii::escape(const std::string& s) {
    std::string ret;
    ret.reserve(s.length()*2);
    for ( std::string::const_iterator it = s.begin(); it != s.end(); ++it ) {
        switch ( *it ) {
        case '\\':
            ret += "\\\\";
            break;
        case '\n':
            ret += "\\|";
            break;
        default:
            ret += *it;
        }
    }
    return ret;
}

void WriterAscii::write_vertex(const ConstGenVertexPtr& v) {
    flush();
    std::string vlist;
    std::vector<int> pids;
    pids.reserve(v->particles_in().size());
    for (const ConstGenParticlePtr& p: v->particles_in()) pids.emplace_back(p->id());
    //We order pids to be able to compare ascii files
    std::sort(pids.begin(), pids.end());
    for (const auto& p: pids) vlist.append( std::to_string(p).append(",") );
    if ( !pids.empty() ) vlist.pop_back();
    const FourVector &pos = v->position();
    if ( !pos.is_zero() ) {
        m_cursor += sprintf(m_cursor, m_vertex_long_printf_specifier.c_str(),  v->id(), v->status(), vlist.c_str(), pos.x(), pos.y(), pos.z(), pos.t() );
    } else {
        m_cursor += sprintf(m_cursor, m_vertex_short_printf_specifier.c_str(), v->id(), v->status(), vlist.c_str());
    }
    flush();
}


inline void WriterAscii::flush() {
    // The maximum size of single add to the buffer (other than by
    // using WriterAscii::write_string) should not be larger than 256. This is a safe value as
    // we will not allow precision larger than 24 anyway
    if ( m_buffer + m_buffer_size < m_cursor + 512 ) {
        std::ptrdiff_t length = m_cursor - m_buffer;
        m_stream->write(m_buffer, length);
        m_cursor = m_buffer;
    }
}


inline void WriterAscii::forced_flush() {
    std::ptrdiff_t length = m_cursor - m_buffer;
    m_stream->write(m_buffer, length);
    m_cursor = m_buffer;
}


void WriterAscii::write_run_info() {
    allocate_buffer();

    // If no run info object set, create a dummy one.
    if ( !run_info() ) set_run_info(std::make_shared<GenRunInfo>());

    const std::vector<std::string> names = run_info()->weight_names();

    if ( !names.empty() ) {
        std::string out = names[0];
        for ( int i = 1, N = names.size(); i < N; ++i ) {
            out += "\n" + names[i];
        }
        m_cursor += sprintf(m_cursor, "W ");
        flush();
        write_string(escape(out));
        m_cursor += sprintf(m_cursor, "\n");
    }

    for (const auto& tool: run_info()->tools()) {
        std::string out = "T " + tool.name + "\n" + tool.version + "\n" + tool.description;
        write_string(escape(out));
        m_cursor += sprintf(m_cursor, "\n");
    }


    for ( const auto& att: run_info()->attributes() ) {
        std::string st;
        if ( !att.second->to_string(st) ) {
            HEPMC3_WARNING_LEVEL(300,"WriterAscii::write_run_info: problem serializing attribute: " << att.first)
        }
        else {
            m_cursor += sprintf(m_cursor, "A ");
            write_string(att.first);
            flush();
            m_cursor += sprintf(m_cursor, " ");
            write_string(escape(st));
            m_cursor += sprintf(m_cursor, "\n");
            flush();
        }
    }
}

void WriterAscii::write_particle(const ConstGenParticlePtr& p, int second_field) {
    flush();
    m_cursor += sprintf(m_cursor, m_particle_printf_specifier.c_str(), p->id(), second_field, p->pid(), p->momentum().px(), p->momentum().py(), p->momentum().pz(), p->momentum().e(), p->generated_mass(), p->status());
    flush();
}


inline void WriterAscii::write_string(const std::string &str) {
    // First let's check if string will fit into the buffer
    if ( m_buffer + m_buffer_size > m_cursor + str.length() ) {
        strncpy(m_cursor, str.data(), str.length());
        m_cursor += str.length();
        flush();
    }
    // If not, flush the buffer and write the string directly
    else {
        forced_flush();
        m_stream->write(str.data(), str.length());
    }
}


void WriterAscii::close() {
    if (!m_stream) return;
    auto* ofs = dynamic_cast<std::ofstream*>(m_stream);
    if (ofs && !ofs->is_open()) return;
    forced_flush();
    const std::string footer("HepMC::Asciiv3-END_EVENT_LISTING\n\n");
    if (m_stream) m_stream->write(footer.data(),footer.length());
    m_stream = nullptr;
    if (ofs) ofs->close();
}
bool WriterAscii::failed() { return (bool)m_file.rdstate(); }

void WriterAscii::set_precision(const int& prec ) {
    if (prec < 2 || prec > 24) return;
    m_precision = prec;
}

int WriterAscii::precision() const {
    return m_precision;
}

void WriterAscii::set_buffer_size(const size_t& size ) {
    if (m_buffer) return;
    if (size < 1024) return;
    m_buffer_size = size;
}


} // namespace HepMC3

// ============================================================
// hepmc3/ReaderAscii.cc
// Copyright (C) 2014-2023 The HepMC collaboration
// ============================================================
// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2023 The HepMC collaboration (see AUTHORS for details)
//
///
/// @file ReaderAscii.cc
/// @brief Implementation of \b class ReaderAscii
///
#include <array>
#include <cstring>
#include <sstream>

#include "HepMC3/ReaderAscii.h"

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Units.h"

namespace HepMC3 {


ReaderAscii::ReaderAscii(const std::string &filename)
    : m_file(filename), m_isstream(false)
{
    if ( !m_file.is_open() ) {
        HEPMC3_ERROR_LEVEL(100,"ReaderAscii: could not open input file: " << filename)
    }
    set_run_info(std::make_shared<GenRunInfo>());
}

ReaderAscii::ReaderAscii(std::istream & stream)
    : m_stream(&stream), m_isstream(true)
{
    if ( !m_stream->good() ) {
        HEPMC3_ERROR_LEVEL(100,"ReaderAscii: could not open input stream ")
    }
    set_run_info(std::make_shared<GenRunInfo>());
}


ReaderAscii::ReaderAscii(std::shared_ptr<std::istream> s_stream)
    : m_shared_stream(s_stream), m_stream(s_stream.get()), m_isstream(true)
{
    if ( !m_stream->good() ) {
        HEPMC3_ERROR_LEVEL(100,"ReaderAscii: could not open input stream ")
    }
    set_run_info(std::make_shared<GenRunInfo>());
}

ReaderAscii::~ReaderAscii() { if (!m_isstream) close(); }

bool ReaderAscii::skip(const int n)
{
    std::array<char, 262144> buf{};
    bool               event_context    = false;
    bool               run_info_context    = false;
    int nn = n;
    while (!failed()) {
        char  peek(0);
        if ( (!m_file.is_open()) && (!m_isstream) ) return false;
        m_isstream ? peek = m_stream->peek() : peek = m_file.peek();
        if ( peek == 'E' ) { event_context = true; nn--; }
        //We have to read each run info.
        if ( !event_context && ( peek == 'W' || peek == 'A' || peek == 'T' ) ) {
            m_isstream ? m_stream->getline(buf.data(), buf.size()) : m_file.getline(buf.data(), buf.size());
            if (!run_info_context) {
                set_run_info(std::make_shared<GenRunInfo>());
                run_info_context = true;
            }
            if ( peek == 'W' ) {
                parse_weight_names(buf.data());
            }
            if ( peek == 'T' ) {
                parse_tool(buf.data());
            }
            if ( peek == 'A' ) {
                parse_run_attribute(buf.data());
            }
        }
        if ( event_context && ( peek == 'V' || peek == 'P' ) ) event_context=false;
        if (nn < 0) return true;
        m_isstream ? m_stream->getline(buf.data(), buf.size()) : m_file.getline(buf.data(), buf.size());
    }
    return true;
}


bool ReaderAscii::read_event(GenEvent &evt) {
    if ( (!m_file.is_open()) && (!m_isstream) ) return false;

    char               peek(0);
    std::array<char, 262144> buf{};
    bool               event_context    = false;
    bool               parsed_weights    = false;
    bool               parsed_particles_or_vertices    = false;
    bool               run_info_context    = false;
    bool               is_parsing_successful  = true;
    std::pair<int, int> vertices_and_particles(0, 0);

    evt.clear();
    evt.set_run_info(run_info());
    m_io_explicit.clear();
    m_io_implicit.clear();
    m_io_implicit_ids.clear();
    m_io_explicit_ids.clear();
    m_data.particles.clear();
    m_data.vertices.clear();
    m_data.links1.clear();
    m_data.links2.clear();
    m_data.attribute_id.clear();
    m_data.attribute_name.clear();
    m_data.attribute_string.clear();
    //
    // Parse event, vertex and particle information
    //
    while (!failed()) {
        m_isstream ? m_stream->getline(buf.data(), buf.size()) : m_file.getline(buf.data(), buf.size());

        if ( strlen(buf.data()) == 0 ) continue;

        // Check for ReaderAscii header/footer
        if ( strncmp(buf.data(), "HepMC", 5) == 0 ) {
            if ( strncmp(buf.data(), "HepMC::Version", 14) != 0 && strncmp(buf.data(), "HepMC::Asciiv3", 14) != 0 )
            {
                HEPMC3_WARNING_LEVEL(500,"ReaderAscii: found unsupported expression in header. Will close the input.")
                std::cout << buf.data() << std::endl;
                m_isstream ? m_stream->clear(std::ios::eofbit) : m_file.clear(std::ios::eofbit);
            }
            if (event_context) {
                is_parsing_successful = true;
                break;
            }
            continue;
        }

        switch (buf[0]) {
        case 'E':
            vertices_and_particles = parse_event_information( buf.data());
            if (vertices_and_particles.second < 0) {
                is_parsing_successful = false;
            } else {
                is_parsing_successful = true;
                event_context   = true;
                parsed_weights = false;
                parsed_particles_or_vertices = false;
            }


            run_info_context   = false;
            break;
        case 'V':
            is_parsing_successful = parse_vertex_information( buf.data());
            parsed_particles_or_vertices =  true;
            break;
        case 'P':
            is_parsing_successful = parse_particle_information( buf.data());
            parsed_particles_or_vertices =  true;
            break;
        case 'W':
            if ( event_context ) {
                is_parsing_successful = parse_weight_values( buf.data());
                parsed_weights=true;
            } else {
                if ( !run_info_context ) {
                    set_run_info(std::make_shared<GenRunInfo>());
                    evt.set_run_info(run_info());
                }
                run_info_context = true;
                is_parsing_successful = parse_weight_names(buf.data());
            }
            break;
        case 'U':
            is_parsing_successful = parse_units( buf.data());
            break;
        case 'T':
            if ( event_context ) {
                //We ignore T in the event context
            } else {
                if ( !run_info_context ) {
                    set_run_info(std::make_shared<GenRunInfo>());
                    evt.set_run_info(run_info());
                }
                run_info_context = true;
                is_parsing_successful = parse_tool(buf.data());
            }
            break;
        case 'A':
            if ( event_context ) {
                is_parsing_successful = parse_attribute( buf.data());
            } else {
                if ( !run_info_context ) {
                    set_run_info(std::make_shared<GenRunInfo>());
                    evt.set_run_info(run_info());
                }
                run_info_context = true;
                is_parsing_successful = parse_run_attribute(buf.data());
            }
            break;
        default:
            HEPMC3_WARNING_LEVEL(500,"ReaderAscii: skipping unrecognised prefix: " << buf[0])
            is_parsing_successful = true;
            break;
        }

        if ( !is_parsing_successful ) break;

        // Check for next event or run info
        m_isstream ? peek = m_stream->peek() : peek = m_file.peek();
        //End of event. The next entry is event.
        if ( event_context &&  peek == 'E' ) break;

        //End of event. The next entry is run info which starts from weight name.
        if ( event_context &&  peek == 'W' && parsed_weights  ) break;

        //End of event. The next entry is run info which starts from attribute.
        if ( event_context &&  peek == 'A' && parsed_particles_or_vertices  ) break;

        //End of event. The next entry is run info which starts from tool.
        if ( event_context &&  peek == 'T' ) break;

    }

    /// Insert the implicit vertices in the gaps of explicit vertices:
    /// Find the gaps looping over the explicit vertices
    int currid = -static_cast<int>(m_data.vertices.size());
    auto fir = m_io_implicit_ids.rbegin();
    for (const auto& iofirst: m_io_explicit_ids) {
        for (; currid < iofirst; ++currid, ++fir) {
            if (fir == m_io_implicit_ids.rend()) {
                HEPMC3_ERROR_LEVEL(600,"ReaderAscii: not enough implicit vertices")
            }
            /// Found a gap in ids, insert an implicit vertex into a list of gaps.
            m_io_explicit[currid] = std::move(m_io_implicit[*fir]);
        }
        ++currid;
    }

    for (const auto& io: m_io_explicit) {
        for (const auto& i: io.second.first) { m_data.links1.push_back(i); m_data.links2.push_back(io.first); }
        for (const auto& o: io.second.second) { m_data.links1.push_back(io.first); m_data.links2.push_back(o); }
    }
    evt.read_data(m_data);

    // Check if all particles and vertices were parsed
    if ((int)evt.particles().size() > vertices_and_particles.second) {
        HEPMC3_ERROR_LEVEL(600,"ReaderAscii: too many particles were parsed")
        printf("%zu  vs  %i expected\n", evt.particles().size(), vertices_and_particles.second);
        is_parsing_successful = false;
    }
    if ((int)evt.particles().size() < vertices_and_particles.second) {
        HEPMC3_ERROR_LEVEL(600,"ReaderAscii: too few  particles were parsed")
        printf("%zu  vs  %i expected\n", evt.particles().size(), vertices_and_particles.second);
        is_parsing_successful = false;
    }

    if ((int)evt.vertices().size()  > vertices_and_particles.first) {
        HEPMC3_ERROR_LEVEL(600,"ReaderAscii: too many vertices were parsed")
        printf("%zu  vs  %i expected\n", evt.vertices().size(), vertices_and_particles.first);
        is_parsing_successful =  false;
    }

    if ((int)evt.vertices().size()  < vertices_and_particles.first) {
        HEPMC3_ERROR_LEVEL(600,"ReaderAscii: too few vertices were parsed")
        printf("%zu  vs  %i expected\n", evt.vertices().size(), vertices_and_particles.first);
        is_parsing_successful =  false;
    }
    // Check if there were HEPMC3_ERRORs during parsing
    if ( !is_parsing_successful ) {
        HEPMC3_ERROR_LEVEL(600,"ReaderAscii: event parsing failed. Returning empty event")
        HEPMC3_DEBUG(1, "Parsing failed at line:" << std::endl << buf.data())

        evt.clear();
        m_isstream ? m_stream->clear(std::ios::badbit) : m_file.clear(std::ios::badbit);

        return false;
    }


    return true;
}


std::pair<int, int> ReaderAscii::parse_event_information(const char *buf) {
    static const std::pair<int, int>  err(-1, -1);
    std::pair<int, int>               ret(-1, -1);
    const char                 *cursor   = buf;
    FourVector&                  position = m_data.event_pos;

    // event number
    if ( !(cursor = strchr(cursor+1, ' ')) ) return err;
    m_data.event_number = atoi(cursor);

    // num_vertices
    if ( !(cursor = strchr(cursor+1, ' ')) ) return err;
    ret.first = atoi(cursor);

    // num_particles
    if ( !(cursor = strchr(cursor+1, ' ')) ) return err;
    ret.second = atoi(cursor);
    m_data.vertices = std::vector<GenVertexData>(ret.first);
    m_data.particles = std::vector<GenParticleData>(ret.second);

    m_data.links1.reserve(ret.second*2);
    m_data.links2.reserve(ret.second*2);
    m_data.attribute_id.reserve(ret.second + ret.first);
    m_data.attribute_name.reserve(ret.second + ret.first);
    m_data.attribute_string.reserve(ret.second + ret.first);
    m_io_implicit_ids.reserve(ret.second);
    // check if there is position information
    if ( (cursor = strchr(cursor+1, '@')) ) {
        // x
        if ( !(cursor = strchr(cursor+1, ' ')) ) return err;
        position.setX(atof(cursor));

        // y
        if ( !(cursor = strchr(cursor+1, ' ')) ) return err;
        position.setY(atof(cursor));

        // z
        if ( !(cursor = strchr(cursor+1, ' ')) ) return err;
        position.setZ(atof(cursor));

        // t
        if ( !(cursor = strchr(cursor+1, ' ')) ) return err;
        position.setT(atof(cursor));
    }

    HEPMC3_DEBUG(10, "ReaderAscii: E: " << m_data.event_number << " (" <<ret.first << "V, " << ret.second << "P)")

    return ret;
}


bool ReaderAscii::parse_weight_values(const char *buf) {
    std::istringstream iss(buf + 1);
    std::vector<double> wts;
    double w = 0.0;
    while (iss >> w) wts.emplace_back(w);
    if ( run_info() && !run_info()->weight_names().empty()
            && run_info()->weight_names().size() != wts.size() ) {
        throw std::logic_error("ReaderAscii::parse_weight_values: "
                               "The number of weights ("+std::to_string((long long int)(wts.size()))+") does not match "
                               "the  number weight names("+std::to_string((long long int)(run_info()->weight_names().size()))+") in the GenRunInfo object");
    }
    m_data.weights = wts;

    return true;
}


bool ReaderAscii::parse_units(const char *buf) {
    const char *cursor = buf;

    // momentum
    if ( !(cursor = strchr(cursor+1, ' ')) ) return false;
    ++cursor;
    m_data.momentum_unit = Units::momentum_unit(cursor);

    // length
    if ( !(cursor = strchr(cursor+1, ' ')) ) return false;
    ++cursor;
    m_data.length_unit = Units::length_unit(cursor);

    HEPMC3_DEBUG(10, "ReaderAscii: U: " << Units::name(m_data.momentum_unit) << " " << Units::name(m_data.length_unit))

    return true;
}


bool ReaderAscii::parse_vertex_information(const char *buf) {
    GenVertexPtr  data = std::make_shared<GenVertex>();
    const char   *cursor          = buf;
    const char   *cursor2         = nullptr;
    int           id              = 0;

    // id
    if ( !(cursor = strchr(cursor+1, ' ')) ) return false;
    id = atoi(cursor);

    // status
    if ( !(cursor = strchr(cursor+1, ' ')) ) return false;
    m_data.vertices[-id-1].status = atoi(cursor);
    FourVector&  position = m_data.vertices[-id-1].position;

    // skip to the list of particles
    if ( !(cursor = strchr(cursor+1, '[')) ) return false;

    while (true) {
        ++cursor;             // skip the '[' or ',' character
        cursor2     = cursor; // save cursor position
        int  particle_in = atoi(cursor);

        // add incoming particle to the vertex
        if (particle_in > 0) {
            //If the particle has not been red yet, we store its id to add the particle later.
            m_io_explicit[id].first.insert(particle_in);
        }

        // check for next particle or end of particle list
        if ( !(cursor = strchr(cursor+1, ',')) ) {
            if ( !(cursor = strchr(cursor2+1, ']')) ) return false;
            break;
        }
    }

    // check if there is position information
    if ( (cursor = strchr(cursor+1, '@')) ) {
        // x
        if ( !(cursor = strchr(cursor+1, ' ')) ) return false;
        position.setX(atof(cursor));

        // y
        if ( !(cursor = strchr(cursor+1, ' ')) ) return false;
        position.setY(atof(cursor));

        // z
        if ( !(cursor = strchr(cursor+1, ' ')) ) return false;
        position.setZ(atof(cursor));

        // t
        if ( !(cursor = strchr(cursor+1, ' ')) ) return false;
        position.setT(atof(cursor));
    }

    return true;
}


bool ReaderAscii::parse_particle_information(const char *buf) {
    const char     *cursor  = buf;
    int             mother_id = 0;

    // verify id
    if ( !(cursor = strchr(cursor+1, ' ')) ) return false;

    int id = atoi(cursor);
    if ( id < 1 || id > static_cast<int>(m_data.particles.size()) ) {
        HEPMC3_ERROR_LEVEL(600,"ReaderAscii: particle ID is out of expected range.")
        return false;
    }

    FourVector&      momentum = m_data.particles[id-1].momentum;
    // mother id
    if ( !(cursor = strchr(cursor+1, ' ')) ) return false;
    mother_id = atoi(cursor);
    if ( mother_id < -static_cast<int>(m_data.vertices.size()) || mother_id > static_cast<int>(m_data.particles.size()) ) {
        HEPMC3_ERROR_LEVEL(600,"ReaderAscii: ID of particle mother is out of expected range.")
        return false;
    }

    if ( mother_id > 0) {
        /// Parent object is a particle, i.e. the vertex is implicit.
        /// If the vertex is not known -- mark its first appearence.
        if (m_io_implicit.count(mother_id) == 0) m_io_implicit_ids.push_back(mother_id);
        m_io_implicit[mother_id].first.insert(mother_id);
        m_io_implicit[mother_id].second.insert(id);
    } else {
        m_io_explicit[mother_id].second.insert(id);
        m_io_explicit_ids.insert(mother_id);
    }
    // pdg id
    if ( !(cursor = strchr(cursor+1, ' ')) ) return false;
    m_data.particles[id-1].pid = atoi(cursor);

    // px
    if ( !(cursor = strchr(cursor+1, ' ')) ) return false;
    momentum.setPx(atof(cursor));

    // py
    if ( !(cursor = strchr(cursor+1, ' ')) ) return false;
    momentum.setPy(atof(cursor));

    // pz
    if ( !(cursor = strchr(cursor+1, ' ')) ) return false;
    momentum.setPz(atof(cursor));

    // pe
    if ( !(cursor = strchr(cursor+1, ' ')) ) return false;
    momentum.setE(atof(cursor));

    // m
    if ( !(cursor = strchr(cursor+1, ' ')) ) return false;
    m_data.particles[id-1].mass = atof(cursor);
    m_data.particles[id-1].is_mass_set = true;

    // status
    if ( !(cursor = strchr(cursor+1, ' ')) ) return false;
    m_data.particles[id-1].status = atoi(cursor);

    return true;
}


bool ReaderAscii::parse_attribute(const char *buf) {
    const char     *cursor  = buf;
    const char     *cursor2 = buf;
    std::array<char, 512> name{};
    int             id = 0;

    if ( !(cursor = strchr(cursor+1, ' ')) ) return false;
    id = atoi(cursor);

    if ( !(cursor  = strchr(cursor+1, ' ')) ) return false;
    ++cursor;

    if ( !(cursor2 = strchr(cursor, ' ')) ) return false;
    snprintf(name.data(), name.size(), "%.*s", (int)(cursor2-cursor), cursor);

    cursor = cursor2+1;

    m_data.attribute_id.push_back(id);
    m_data.attribute_name.emplace_back(name.data());
    m_data.attribute_string.push_back(unescape(cursor));

    return true;
}

bool ReaderAscii::parse_run_attribute(const char *buf) {
    const char     *cursor  = buf;
    const char     *cursor2 = buf;
    std::array<char, 512> name{};

    if ( !(cursor = strchr(cursor+1, ' ')) ) return false;
    ++cursor;

    if ( !(cursor2 = strchr(cursor, ' ')) ) return false;
    snprintf(name.data(), name.size(), "%.*s", (int)(cursor2-cursor), cursor);

    cursor = cursor2+1;

    std::shared_ptr<StringAttribute> att =
        std::make_shared<StringAttribute>(StringAttribute(unescape(cursor)));

    run_info()->add_attribute(std::string(name.data()), att);

    return true;
}


bool ReaderAscii::parse_weight_names(const char *buf) {
    const char     *cursor  = buf;

    if ( !(cursor = strchr(cursor+1, ' ')) ) return false;
    ++cursor;

    std::istringstream iss(unescape(cursor));
    std::vector<std::string> names;
    std::string name;
    while (iss >> name) names.emplace_back(name);

    run_info()->set_weight_names(names);

    return true;
}

bool ReaderAscii::parse_tool(const char *buf) {
    const char     *cursor  = buf;

    if ( !(cursor = strchr(cursor+1, ' ')) ) return false;
    ++cursor;
    std::string line = unescape(cursor);
    GenRunInfo::ToolInfo tool;
    std::string::size_type pos = line.find('\n');
    tool.name = line.substr(0, pos);
    line = line.substr(pos + 1);
    pos = line.find('\n');
    tool.version = line.substr(0, pos);
    tool.description = line.substr(pos + 1);
    run_info()->tools().emplace_back(tool);

    return true;
}


std::string ReaderAscii::unescape(const std::string& s) {
    std::string ret;
    ret.reserve(s.length());
    for ( std::string::const_iterator it = s.begin(); it != s.end(); ++it ) {
        if ( *it == '\\' ) {
            ++it;
            if ( *it == '|' ) {
                ret += '\n';
            }
            else {
                ret += *it;
            }
        } else
        {ret += *it;}
    }

    return ret;
}

bool ReaderAscii::failed() { return m_isstream ? (bool)m_stream->rdstate() :(bool)m_file.rdstate(); }

void ReaderAscii::close() {
    if ( !m_file.is_open()) return;
    m_file.close();
}


} // namespace HepMC3

#endif /* MARLEY_FOUND_HEPMC3 */
