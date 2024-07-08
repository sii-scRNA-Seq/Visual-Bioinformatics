import { By } from '@angular/platform-browser';
import { ComponentFixture, TestBed } from '@angular/core/testing';
import { HttpClientTestingModule } from '@angular/common/http/testing';
import { MatCardModule } from '@angular/material/card';
import { MatProgressSpinnerModule } from '@angular/material/progress-spinner';
import { MatSnackBarModule } from '@angular/material/snack-bar';

import { BlockService } from '../block.service';
import { CanvasComponent } from './canvas.component';
import { CodeBlockComponent } from '../code-block/code-block.component';
import { MockBlockService } from '../mock-block.service';
import { OutputService } from '../output.service';
import { MockOutputService } from '../mock-output.service';

describe('CanvasComponent', () => {
  let component: CanvasComponent;
  let fixture: ComponentFixture<CanvasComponent>;

  beforeEach(() => {
    TestBed.configureTestingModule({
      declarations: [
        CanvasComponent,
        CodeBlockComponent
      ],
      imports: [
        HttpClientTestingModule,
        MatCardModule,
        MatProgressSpinnerModule,
        MatSnackBarModule,
      ],
      providers: [
        { provide: BlockService, useClass: MockBlockService },
        { provide: OutputService, useClass: MockOutputService },
      ],
    });
    fixture = TestBed.createComponent(CanvasComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });

  describe('Running blocks', () => {
    it ('should call blockService.executeBlocks when run button is clicked', () => {
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'executeBlocks');
      const button = fixture.debugElement.query(By.css('button'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.executeBlocks).toHaveBeenCalledTimes(1);
    });

    it('should be disabled when no blocks have been added', () => {
      component.executingBlocks = false;
      component.blockList = [];
      fixture.detectChanges();
      expect(fixture.debugElement.query(By.css('button')).nativeElement.disabled).toEqual(true);
    });

    it('should not be disabled when blocks have been added', () => {
      component.executingBlocks = false;
      component.blockList = [{
        blockId: 'loaddata',
        blockUUID: '',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [],
      }];
      fixture.detectChanges();
      expect(fixture.debugElement.query(By.css('button')).nativeElement.disabled).toEqual(false);
    });

    it('should be disabled when all blocks have been removed', () => {
      component.executingBlocks = false;
      component.blockList = [];
      fixture.detectChanges();
      component.blockList = [{
        blockId: 'loaddata',
        blockUUID: '',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [],
      }];
      fixture.detectChanges();
      component.blockList = [];
      fixture.detectChanges();
      expect(fixture.debugElement.query(By.css('button')).nativeElement.disabled).toEqual(true);
    });

    it('should be available when blocks are not being executed', () => {
      component.executingBlocks = false;
      fixture.detectChanges();
      expect(fixture.debugElement.query(By.css('button'))).not.toBeNull();
    });
    
    it('should not be available while blocks are being executed', () => {
      component.executingBlocks = false;
      fixture.detectChanges();
      expect(fixture.debugElement.query(By.css('button'))).not.toBeNull();
      component.executingBlocks = true;
      fixture.detectChanges();
      expect(fixture.debugElement.query(By.css('button'))).toBeNull();
    });

    it ('should become available once blocks have stopped being executed', () => {
      component.executingBlocks = true;
      fixture.detectChanges();
      expect(fixture.debugElement.query(By.css('button'))).toBeNull();
      component.executingBlocks = false;
      fixture.detectChanges();
      expect(fixture.debugElement.query(By.css('button'))).not.toBeNull();
    });
  });

  describe('Progress spinner', () => {
    it('should not appear when blocks are not being executed', () => {
      component.executingBlocks = false;
      fixture.detectChanges();
      expect(fixture.debugElement.query(By.css('mat-spinner'))).toBeNull();
    });

    it('should appear while blocks are being executed', () => {
      component.executingBlocks = false;
      fixture.detectChanges();
      expect(fixture.debugElement.query(By.css('mat-spinner'))).toBeNull();
      component.executingBlocks = true;
      fixture.detectChanges(); 
      expect(fixture.debugElement.query(By.css('mat-spinner'))).not.toBeNull();
    });

    it('should disappear once blocks have stopped being executed', () => {
      component.executingBlocks = true;
      fixture.detectChanges(); 
      expect(fixture.debugElement.query(By.css('mat-spinner'))).not.toBeNull();
      component.executingBlocks = false;
      fixture.detectChanges();
      expect(fixture.debugElement.query(By.css('mat-spinner'))).toBeNull();
    });
  });
});
