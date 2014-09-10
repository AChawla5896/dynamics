#ifndef XML_TEST_H
#define XML_TEST_H

#include <util/misc/XmlBase.h>
#include <util/misc/XmlAttribute.h>
#include <util/misc/XmlStartTag.h>
#include <util/misc/XmlEndTag.h>
#include <util/global.h>

#ifdef UTIL_MPI
#ifndef TEST_MPI
#define TEST_MPI
#endif
#endif

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;

class XmlTest : public UnitTest 
{

public:

   void setUp()
   {};

   void tearDown()
   {};

   void testXmlBaseSetString() 
   {
      printMethod(TEST_FUNC);

      XmlBase parser;
      std::string string("this is the string");
      parser.setString(string, 0);
      TEST_ASSERT(parser.string() == string);
      TEST_ASSERT(parser.cursor() == 0);
      TEST_ASSERT(parser.c() == 't');
   }

   void testXmlBaseMove() 
   {
      printMethod(TEST_FUNC);

      XmlBase parser;
      std::string string("   this is the string");
      parser.setString(string, 0);
      TEST_ASSERT(parser.string() == string);
      TEST_ASSERT(parser.cursor() == 0);
      TEST_ASSERT(parser.c() == ' ');
      parser.skip();
      TEST_ASSERT(parser.cursor() == 3);
      TEST_ASSERT(parser.c() == 't');
      parser.next();
      TEST_ASSERT(parser.cursor() == 4);
      TEST_ASSERT(parser.c() == 'h');
      while (!parser.isEnd()) {
         parser.next();
      }
      TEST_ASSERT(parser.isEnd());
      TEST_ASSERT(parser.c() == '\0');
      TEST_ASSERT(parser.cursor() == 21);
      
   }

   void testXmlAttributeMatch1() 
   {
      printMethod(TEST_FUNC);

      XmlAttribute parser;
      std::string string("  this =\"35.5\"");
      bool result = parser.match(string, 0);
      TEST_ASSERT(result == true);
      //std::cout << "\n" << parser.label();
      //std::cout << "\n" << parser.value().str();
      double value;
      parser.value() >>  value;
      TEST_ASSERT(parser.label() == std::string("this"));
      TEST_ASSERT(value == 35.5);
   }

   void testXmlAttributeMatch2() 
   {
      printMethod(TEST_FUNC);

      XmlBase parent;
      XmlAttribute parser;
      std::string string("  this =\"35.5\"u");
      parent.setString(string, 0);
      bool result = parser.match(parent);
      TEST_ASSERT(result == true);
      //std::cout << "\n" << parser.label();
      //std::cout << "\n" << parser.value().str();
      double value;
      parser.value() >>  value;
      TEST_ASSERT(parser.label() == std::string("this"));
      TEST_ASSERT(value == 35.5);
      TEST_ASSERT(parser.cursor() == parent.cursor());
      TEST_ASSERT(parent.c() == 'u');
      TEST_ASSERT(parent.cursor() == 14);
   }

   void testXmlStartTag() 
   {
      printMethod(TEST_FUNC);

      XmlStartTag tag;
      std::string string("<Label double =\"35.5\" string= \"glib\" />");
      bool result = tag.matchLabel(string, 0);
      TEST_ASSERT(result);
      TEST_ASSERT(tag.label() == std::string("Label"));
      std::cout << std::endl << "label = " << tag.label();
      XmlAttribute attribute;
      while (tag.matchAttribute(attribute)) {
         std::cout << std::endl 
                   << " label = " << attribute.label() 
                   << " value = " << attribute.value().str();
         if (attribute.label() == "double") {
            double value;
            attribute.value() >> value;
            //std::cout << std::endl << "Extracted value";
            TEST_ASSERT(value == 35.5);
         }
      }
      std::cout << std::endl;
      TEST_ASSERT(tag.endBracket());
   }

   void testXmlEndTag() 
   {
      printMethod(TEST_FUNC);
      XmlEndTag tag;
      std::string string("</Thing>");
      bool result = tag.match(string, 0);
      TEST_ASSERT(result);
      TEST_ASSERT(tag.label() == "Thing");

      string = "</ Knobby >";
      result = tag.match(string, 0);
      TEST_ASSERT(result);
      TEST_ASSERT(tag.label() == "Knobby");
   }
};

TEST_BEGIN(XmlTest)
TEST_ADD(XmlTest, testXmlBaseSetString)
TEST_ADD(XmlTest, testXmlBaseMove)
TEST_ADD(XmlTest, testXmlAttributeMatch1)
TEST_ADD(XmlTest, testXmlAttributeMatch2)
TEST_ADD(XmlTest, testXmlStartTag)
TEST_ADD(XmlTest, testXmlEndTag)
TEST_END(XmlTest)

#endif