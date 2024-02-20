import { ComponentFixture, TestBed } from '@angular/core/testing';
import { DomSanitizer } from '@angular/platform-browser';
import { HttpClientTestingModule } from '@angular/common/http/testing';
import { MatCardModule } from '@angular/material/card';
import { MatSnackBarModule } from '@angular/material/snack-bar';

import { Block } from '../block.interface';
import { MockOutputService } from '../mock-output.service';
import { OutputDisplayComponent } from './output-display.component';
import { OutputService } from '../output.service';

describe('OutputDisplayComponent', () => {
  let component: OutputDisplayComponent;
  let fixture: ComponentFixture<OutputDisplayComponent>;
  let sanitizer: DomSanitizer;


  beforeEach(() => {
    TestBed.configureTestingModule({
      declarations: [OutputDisplayComponent],
      imports: [
        HttpClientTestingModule,
        MatCardModule,
        MatSnackBarModule,
      ],
      providers: [
        { provide: OutputService, useClass: MockOutputService }
      ],
    });
    fixture = TestBed.createComponent(OutputDisplayComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
    sanitizer = TestBed.inject(DomSanitizer);
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });

  describe('OutputList', () => {
    it('should have text outputs remain the same', () => {
      const block: Block = {
        blockId: 'loaddata',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [],
      };
      const outputService: OutputService = TestBed.inject(OutputService);
      outputService.executeBlock(block);
      expect(component.outputList[0].text == 'Some text here');
    });

    it('should replace image outputs with a sanitised SafeUrl', () => {
      const block: Block = {
        blockId: 'loaddata',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [],
      };
      const outputService: OutputService = TestBed.inject(OutputService);
      const spy = spyOn(sanitizer, 'bypassSecurityTrustUrl');
      outputService.executeBlock(block);
      expect(spy).toHaveBeenCalledTimes(1);
      const expectedValue = sanitizer.bypassSecurityTrustUrl('data:image/png;base64,' + 'An image here');
      expect(component.outputList[1].img).toBe(expectedValue);
    });
  });
});
