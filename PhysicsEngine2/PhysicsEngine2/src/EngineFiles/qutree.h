#ifndef QUTREE_H
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <vector>
#include <array>
#include <list>
#include <memory>
#include <functional>
#include "Vector2D.h"
#include "QTRstruct.h"
const int MAX_QDEPTH = 8;


//the javidx9 tutorial
template <typename T>
struct QuadTreeItemLocation
{
    typename std::list<std::pair<QTRect, T>>* container = nullptr;
    typename std::list<std::pair<QTRect, T>>::iterator iterator;
};

template <typename pT>
class DynamicQuadTree
{
public:
    DynamicQuadTree(const QTRect& size = {Vector2D(-1,-1),Vector2D(2,2)}, const size_t nDepth = 0, const size_t nMaxDepth = 6)
    {
        m_depth = nDepth;
        m_rect = size;
        m_maxdepth = nMaxDepth;
        resize(m_rect);
    }

    void drawQT(const QTRect& rect) {
        float x1 = rect.pos.x;
        float y1 = rect.pos.y;
        float x2 = rect.pos.x + rect.dimensions.x;
        float y2 = rect.pos.y + rect.dimensions.y;

        glBegin(GL_LINE_LOOP);
        glColor3f(1.0f, 1.0f, 1.0f);
        glVertex2f(x1, y1);
        glVertex2f(x2, y1);
        glVertex2f(x2, y2);
        glVertex2f(x1, y2);
        glEnd();
    }

    // Insert a region into this area
    QuadTreeItemLocation<pT> insert(const pT item, const QTRect& itemsize)
    {
        for (int i = 0; i < 4; i++)
        {
            if (m_rChild[i].containr(itemsize))
            {
                // Have we reached depth limit?
                if (m_depth + 1 < m_maxdepth)
                {
                    // No, so does child exist?
                    if (!m_pChild[i])
                    {
                        // No, so create it
                        m_pChild[i] = std::make_shared<DynamicQuadTree<pT>>(m_rChild[i], m_depth + 1, m_maxdepth);
                    }

                    // Yes, so add item to it
                    return m_pChild[i]->insert(item, itemsize);
                }
            }
        }

        // It didnt fit, so item must belong to this geom2d::rect<CTYPE>
        m_pItems.push_back({ itemsize, item });
        return { &m_pItems, std::prev(m_pItems.end()) };
    }

    void relocate(pT item, const QTRect& rArea)
    {
        // Remove it
        remove(item);

        // Reinsert it with new location
        insert(item, rArea);
    }

    size_t size() const
    {
        size_t nCount = m_pItems.size();
        for (int i = 0; i < 4; i++)
            if (m_pChild[i]) nCount += m_pChild[i]->size();
        return nCount;
    }
    void searchC(const QTCircle& rArea, std::list<pT>& listItems) const {

    }

    void search(const QTRect& rArea, std::list<pT>& listItems) const
    {
        // First, check for items belonging to this area, add them to the list
        // if there is overlap
        for (const auto& p : m_pItems)
        {
            if (rArea.overlaps(p.first))
                listItems.push_back(p.second);
        }

        // Second, recurse through children and see if they can
        // add to the list
        for (int i = 0; i < 4; i++)
        {
            if (m_pChild[i])
            {
                // If child is entirely contained within area, recursively
                // add all of its children, no need to check boundaries
                if (rArea.containr(m_rChild[i]))
                    m_pChild[i]->items(listItems);

                // If child overlaps with search area then checks need
                // to be made
                else if (m_rChild[i].overlaps(rArea))
                    m_pChild[i]->search(rArea, listItems);
            }
        }


    }

    void items(std::list<pT>& listItems) const
    {
        // No questions asked, just return child items
        for (const auto& p : m_pItems)
            listItems.push_back(p.second);

        for (int i = 0; i < 4; i++)
            if (m_pChild[i]) m_pChild[i]->items(listItems);
    }

    void clear()
    {
        m_pItems.clear();

        for (int i = 0; i < 4; i++)
        {
            if (m_pChild[i])
                m_pChild[i]->clear();
            else
                m_pChild[i].reset();
        }
    }

    void resize(const QTRect& rsSize) {
        clear();
        m_rect = rsSize;
        Vector2D qtChildsize = m_rect.dimensions / 2.0f;
        QTRect q1;
        q1.pos = m_rect.pos;
        q1.dimensions = qtChildsize;
        QTRect q2;
        q2.pos.x = m_rect.pos.x + qtChildsize.x;
        q2.pos.y = m_rect.pos.y;
        q2.dimensions = qtChildsize;
        QTRect q3;
        q3.pos.x = m_rect.pos.x;
        q3.pos.y = m_rect.pos.y + qtChildsize.y;
        q3.dimensions = qtChildsize;
        QTRect q4;
        q4.pos = m_rect.pos + qtChildsize;
        q4.dimensions = qtChildsize;
        m_rChild = { q1,q2,q3,q4 };
    }


    void visualizeQuadTree()
    {

        // Draw the current quadtree node (rectangle)
        drawQT(m_rect);
        //drawQT();
        // Draw the items contained in the current node

        // Recursively draw the children
        for (int i = 0; i < 4; i++)
        {
            if (m_pChild[i])
            {
                m_pChild[i]->visualizeQuadTree();
            }
        }
    }


    const QTRect& area()
    {
        return m_rect;
    }

protected:
    size_t m_depth = 0;
    size_t m_maxdepth = 8;

    // Area of this quadnode
    QTRect m_rect;

    // 4 child areas of this quadnode
    std::array<QTRect, 4> m_rChild{};

    // 4 potential children of this quadnode
    std::array<std::shared_ptr<DynamicQuadTree<pT>>, 4> m_pChild{};

    // Items which belong to this quadnode
    std::list<std::pair<QTRect, pT>> m_pItems;

};

template<typename T>
struct QuadTreeItem
{
    // The item Itself
    T item;
    QTRect boundingbox;
    // A "locator" to the container/iterator that points to this item's iterator in the
    // top level list - phew
    QuadTreeItemLocation<typename std::list<QuadTreeItem<T>>::iterator> pItem;
};

template <typename T>
class QuadTreeContainer
{
    using IQuadtreeContainer = std::list<QuadTreeItem<T>>;

public:
    QuadTreeContainer(const QTRect& size = {Vector2D(-1,-1),Vector2D(2,2)}, const size_t nDepth = 0, const size_t nMaxDepth = 8) : root(size, nDepth, nMaxDepth)
    {

    }

    // Sets the spatial coverage area of teh quadtree
    void resize(const QTRect& rArea)
    {
        root.resize(rArea);
    }

    // Inserts an item into the quadtree
    void insert(const T& item, const QTRect& itemsize)
    {
        QuadTreeItem<T> newItem;
        newItem.item = item;
        newItem.boundingbox = itemsize;
        // Item i stored in container
        m_allItems.emplace_back(newItem);

        // Pointer/Area of item is stored in geom2d::rect<CTYPE> tree
        m_allItems.back().pItem = root.insert(std::prev(m_allItems.end()), itemsize);
    }

    // Returns a std::list of pointers to items within the search area
    std::list<typename IQuadtreeContainer::iterator> search(const QTRect& rArea) const
    {
        std::list<typename IQuadtreeContainer::iterator> listItemPointers;
        root.search(rArea, listItemPointers);
        return listItemPointers;
    }

    void remove(typename IQuadtreeContainer::iterator& item)
    {
        // Iterator points to a QuadTreeItem
        item->pItem.container->erase(item->pItem.iterator);

        // Remove actual item from container
        m_allItems.erase(item);
    }

    void relocate(typename IQuadtreeContainer::iterator& item, const QTRect& itemsize)
    {
        // Remove pointer to item from whichever container its stored in
        QuadTreeItem<T>& qti = *item;

        if (qti.boundingbox != itemsize) {
            item->pItem.container->erase(item->pItem.iterator);

            // Update the items pointer by reinsertion into geom2d::rect<CTYPE> tree
            item->pItem = root.insert(item, itemsize);
        }
        else {
            return;
        }

    }

    typename IQuadtreeContainer::iterator begin()
    {
        return m_allItems.begin();
    }

    typename IQuadtreeContainer::iterator end()
    {
        return m_allItems.end();
    }

    typename IQuadtreeContainer::const_iterator cbegin()
    {
        return m_allItems.cbegin();
    }

    typename IQuadtreeContainer::const_iterator cend()
    {
        return m_allItems.cend();
    }

    size_t size() const
    {
        return root.size();
    }

    void clear()
    {
        m_allItems.clear();
        root.clear();
    }

    const QTRect& area()
    {
        return root.area();
    }

    void visualizeQuadTree()
    {

        return root.visualizeQuadTree();
    }

protected:
    DynamicQuadTree<typename IQuadtreeContainer::iterator> root;
    IQuadtreeContainer m_allItems;
};

#endif //QUTREE_H!


